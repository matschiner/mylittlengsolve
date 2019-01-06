
/*
  
  Solver for the linear hyperbolic equation

  du/dt  +  div (b u) = 0

  by an explicit time-stepping method

*/



#include <solve.hpp>
#include <python_ngstd.hpp>


using namespace ngsolve;
using ngfem::ELEMENT_TYPE;


template<int D>
class Convection {
protected:
    shared_ptr<L2HighOrderFESpace> fes;
    shared_ptr<CoefficientFunction> cfflow;
    FlatVector<int> distant_procs;

    class FacetData {
    public:
        size_t elnr[2];
        int facetnr[2];
        Vector<> flown;

        FacetData() = default;

        FacetData(size_t nip) : flown(nip) { ; }
    };


    class ElementData {
    public:
        MatrixFixWidth<D> flowip;

        ElementData() = default;

        ElementData(size_t ndof, size_t nip) : flowip(nip) { ; }
    };

    Array<FacetData> facetdata;         // new move semantics
    Array<ElementData> elementdata;     // new move semantics

public:

    Convection(shared_ptr<FESpace> afes, shared_ptr<CoefficientFunction> aflow)
            : fes(dynamic_pointer_cast<L2HighOrderFESpace>(afes)), cfflow(aflow) {
        LocalHeap lh(1000000);
        shared_ptr<MeshAccess> ma = fes->GetMeshAccess();
        elementdata.SetAllocSize(ma->GetNE());


        cout << "my id is: " << MyMPI_GetId() << "\nrunning with: " << MyMPI_GetNTasks() << " tasks" << endl;

        if (!fes->AllDofsTogether())
            throw Exception("mlngs-Convection needs 'all_dofs_together=True' for L2-FESpace");

        for (auto ei : ma->Elements()) {
            HeapReset hr(lh);

            auto &fel = dynamic_cast<const DGFiniteElement<D> &> (fes->GetFE(ei, lh));
            const IntegrationRule ir(fel.ElementType(), 2 * fel.Order());

            const_cast<DGFiniteElement<D> &> (fel).PrecomputeShapes(ir);
            const_cast<DGFiniteElement<D> &> (fel).PrecomputeTrace();

            auto &trafo = ma->GetTrafo(ei, lh);
            MappedIntegrationRule<D, D> mir(ir, trafo, lh);

            ElementData edi(fel.GetNDof(), ir.Size());
            cfflow->Evaluate(mir, edi.flowip);
            for (size_t j = 0; j < ir.Size(); j++) {
                Vec<D> flow = mir[j].GetJacobianInverse() * edi.flowip.Row(j);
                flow *= mir[j].GetWeight();     // weight times Jacobian
                edi.flowip.Row(j) = flow;
            }

            elementdata.Append(move(edi));
        }


        Array<int> elnums, fnums, vnums;

        facetdata.SetAllocSize(ma->GetNFacets());

        for (auto i : Range(ma->GetNFacets())) {
            HeapReset hr(lh);

            const DGFiniteElement<D - 1> &felfacet = dynamic_cast<const DGFiniteElement<D - 1> &> (fes->GetFacetFE(i, lh));
            IntegrationRule ir(felfacet.ElementType(), 2 * felfacet.Order());
            const_cast<DGFiniteElement<D - 1> &> (felfacet).PrecomputeShapes(ir);

            FlatVector<int> distant_procs(ma->GetNE(), lh);

            FacetData fai(ir.Size());

            ma->GetFacetElements(i, elnums);

            fai.elnr[1] = -1;
            for (size_t j : Range(elnums)) {
                fai.elnr[j] = elnums[j];
                auto fnums = ma->GetElFacets(ElementId(VOL, elnums[j]));
                fai.facetnr[j] = fnums.Pos(i);
            }


            ELEMENT_TYPE eltype = ma->GetElType(ElementId(VOL, elnums[0]));

            vnums = ma->GetElVertices(ElementId(VOL, elnums[0]));
            Facet2ElementTrafo transform(eltype, vnums);
            FlatVec<D> normal_ref = ElementTopology::GetNormals(eltype)[fai.facetnr[0]];

            size_t nip = ir.Size();

            // transform facet coordinates to element coordinates
            IntegrationRule &irt = transform(fai.facetnr[0], ir, lh);
            MappedIntegrationRule<D, D> mir(irt, ma->GetTrafo(ElementId(VOL, elnums[0]), lh), lh);

            FlatMatrixFixWidth<D> flowir(nip, lh);
            cfflow->Evaluate(mir, flowir);

            for (size_t j = 0; j < nip; j++) {
                Vec<D> normal = Trans(mir[j].GetJacobianInverse()) * normal_ref;

                fai.flown(j) = InnerProduct(normal, flowir.Row(j));
                fai.flown(j) *= ir[j].Weight() * mir[j].GetJacobiDet();
            }

            facetdata.Append(move(fai));
        }

        //        cout << "end on "<<MyMPI_GetId()<<endl;
    }


    void Apply(BaseVector &_vecu, BaseVector &_conv) {
        static Timer t("Convection::Apply");
        //cout << "starting apply" << endl;
        RegionTimer reg(t);
        LocalHeap lh(1000 * 1000);
        LocalHeap gh(1000 * 1000);

        auto vecu = _vecu.FV<double>();
        auto conv = _conv.FV<double>();

        auto ma = fes->GetMeshAccess();
        //pardofs=fes->GetParallelDofs();


        //find out the exchange dofs
        int ranksize = MyMPI_GetNTasks();
        int rankid = MyMPI_GetId();

        Vector<int> shared_facets_counts(ranksize);
        shared_facets_counts = 0;

        for (auto i: Range(ma->GetNFacets())) {
            auto dp = ma->GetDistantProcs(NodeId(StdNodeType(NT_FACET, ma->GetDimension()), i));
            if (dp.Size() == 1) {
                shared_facets_counts[dp[0]]++;
            }
        }

        Array<Vector<double>> data_get[ranksize];
        FlatArray < FlatArray < double >> data_send[ranksize];
        FlatArray < FlatArray < double >> data_get2[ranksize];
        for (auto i: Range(ranksize)) {
            data_get[i].SetAllocSize(shared_facets_counts[i]);
            data_send[i].Assign(shared_facets_counts[i], gh);

            data_get2[i].Assign(shared_facets_counts[i], gh);
        }

        ParallelFor
                (Range(ma->GetNE()), [&](size_t i) {
                    LocalHeap slh = lh.Split(), &lh = slh;

                    auto &fel = static_cast<const ScalarFiniteElement<D> &> (fes->GetFE(ElementId(VOL, i), lh));
                    const IntegrationRule ir(fel.ElementType(), 2 * fel.Order());


                    FlatMatrixFixWidth<D> flowip = elementdata[i].flowip;

                    /*
                     // use this for time-dependent flow (updated version not yet tested)
                     MappedIntegrationRule<D,D> mir(ir, ma->GetTrafo (i, 0, lh), lh);
                     FlatMatrixFixWidth<D> flowip(mir.Size(), lh);
                     cfflow -> Evaluate (mir, flowip);
                     for (size_t j = 0; j < ir.Size(); j++)
                     {
                     Vec<D> flow = mir[j].GetJacobianInverse() * flowip.Row(j);
                     flow *= mir[j].GetWeight();
                     flowip.Row(j) = flow;
                     }
                     */

                    IntRange dn = fes->GetElementDofs(i);

                    size_t nipt = ir.Size();
                    FlatVector<> elui(nipt, lh);
                    FlatMatrixFixWidth<D> flowui(nipt, lh);

                    fel.Evaluate(ir, vecu.Range(dn), elui);

                    for (auto k : Range(nipt))
                        flowui.Row(k) = elui(k) * flowip.Row(k);

                    fel.EvaluateGradTrans(ir, flowui, conv.Range(dn));
                });

        static mutex add_mutex;

        MyMPI_Barrier();
        if (rankid == 0) cout << "finished element data" << endl;


        Vector<int> k(ranksize);
        k = 0;

        ParallelFor
                (Range(ma->GetNFacets()), [&](size_t i) {

                    //        for (auto i: Range(ma->GetNFacets())) {
                    LocalHeap slh = lh.Split(), &lh = slh;

                    const FacetData &fai = facetdata[i];


                    if (fai.elnr[1] != -1) {
                        // internal facet
                        const DGFiniteElement<D> &fel1 = static_cast<const DGFiniteElement<D> &> (fes->GetFE(ElementId(VOL, fai.elnr[0]), lh));
                        const DGFiniteElement<D> &fel2 = static_cast<const DGFiniteElement<D> &> (fes->GetFE(ElementId(VOL, fai.elnr[1]), lh));
                        const DGFiniteElement<D - 1> &felfacet = static_cast<const DGFiniteElement<D - 1> & > (fes->GetFacetFE(i, lh));

                        IntRange dn1 = fes->GetElementDofs(fai.elnr[0]);
                        IntRange dn2 = fes->GetElementDofs(fai.elnr[1]);

                        size_t ndoffacet = felfacet.GetNDof();
                        size_t ndof1 = fel1.GetNDof();
                        size_t ndof2 = fel2.GetNDof();

                        FlatVector<> aelu1(ndof1, lh), aelu2(ndof2, lh);
                        FlatVector<> trace1(ndoffacet, lh), trace2(ndoffacet, lh);

                        fel1.GetTrace(fai.facetnr[0], vecu.Range(dn1), trace1);
                        fel2.GetTrace(fai.facetnr[1], vecu.Range(dn2), trace2);

                        IntegrationRule ir(felfacet.ElementType(), 2 * felfacet.Order());
                        size_t nip = ir.Size();

                        FlatVector<> flown = fai.flown;

                        FlatVector<> tracei1(nip, lh), tracei2(nip, lh);
                        FlatVector<> tracei(nip, lh);

                        felfacet.Evaluate(ir, trace1, tracei1);
                        felfacet.Evaluate(ir, trace2, tracei2);

                        for (size_t j = 0; j < nip; j++)
                            tracei(j) = flown(j) * ((flown(j) > 0) ? tracei1(j) : tracei2(j));

                        felfacet.EvaluateTrans(ir, tracei, trace1);
                        fel1.GetTraceTrans(fai.facetnr[0], trace1, aelu1);
                        fel2.GetTraceTrans(fai.facetnr[1], trace1, aelu2);

                        {
                            lock_guard<mutex> guard(add_mutex);
                            conv.Range(dn1) -= aelu1;
                            conv.Range(dn2) += aelu2;
                        }
                    } else {


                        // boundary facet
                        const DGFiniteElement<D> &fel1 = dynamic_cast<const DGFiniteElement<D> &> (fes->GetFE(ElementId(VOL, fai.elnr[0]), lh));
                        const DGFiniteElement<D - 1> &felfacet = dynamic_cast<const DGFiniteElement<D - 1> & > (fes->GetFacetFE(i, lh));

                        IntRange dn1 = fes->GetElementDofs(fai.elnr[0]);

                        size_t ndoffacet = felfacet.GetNDof();
                        size_t ndof1 = fel1.GetNDof();

                        FlatVector<> elu1(ndof1, lh);
                        FlatVector<> trace1(ndoffacet, lh);

                        fel1.GetTrace(fai.facetnr[0], vecu.Range(dn1), trace1);

                        IntegrationRule ir(felfacet.ElementType(), 2 * felfacet.Order());
                        size_t nip = ir.Size();

                        FlatVector<> flown = fai.flown;
                        FlatVector<> tracei1(nip, lh), tracei(nip, lh);

                        felfacet.Evaluate(ir, trace1, tracei1);

                        auto dp = ma->GetDistantProcs(NodeId(StdNodeType(NT_FACET, ma->GetDimension()), i));
                        if (dp.Size() == 1) {
                            FlatArray < double > tmp(nip, gh);
                            tmp = tracei1;

                            data_send[dp[0]][k[dp[0]]++].Assign(tmp);
                            data_get2[dp[0]][k[dp[0]]++].Assign(tmp);
                        } else {
                            for (size_t j = 0; j < nip; j++)
                                tracei(j) = flown(j) * ((flown(j) > 0) ? tracei1(j) : 0);

                            felfacet.EvaluateTrans(ir, tracei, trace1);
                            fel1.GetTraceTrans(fai.facetnr[0], trace1, elu1);

                            {
                                lock_guard<mutex> guard(add_mutex);
                                conv.Range(dn1) -= elu1;
                            }
                        }
                    }
                });
        /*
        if (rankid==0) cout << "start sending"<<endl;
        Array<MPI_Request> requests;
        for (auto i : Range(ranksize)) {
            if (i == rankid) {
                continue;
            }
            for (auto j : Range(shared_facets_counts[i])) {

                auto t1=MyMPI_ISend(data_send[i][j], i, j);
                requests.Append(t1);

                auto t2 = MyMPI_IRecv(data_get2[i][j], i, j);
                requests.Append(t2);
            }
        }
        for (auto i : Range(ranksize)) {
            cout << data_get2[i]<<endl;
        }

        if (requests.Size() != 0) {
//            cout << "waiting for " << requests.Size() << endl;
            MyMPI_WaitAll(requests);
        }
        MyMPI_Barrier();
        */



        if (rankid == 0) {
            cout << "finish" << endl;
        }

        for (auto i : Range(ranksize)) {
            for (auto j: Range(ranksize)) {

                if (i == j) {
                    continue;
                }
                //                cout << "sending from "<<i<< " to "<<j<<endl;
                if (rankid == i) {
                    for (auto k : Range(shared_facets_counts[j])) {
                        MyMPI_Send(data_send[j][k], j);
                    }
                }

                if (rankid == j) {
                    for (auto k : Range(shared_facets_counts[i])) {
                        Vector<double> tmp1(6);
                        FlatArray < double > tmp2(6, gh);
                        MyMPI_Recv(tmp2, i);
                        data_get2[i][k].Assign(tmp2);
                    }
                }
            }
        }
        for (auto i : Range(ranksize)) {
            cout << data_get2[i] << endl;
        }

        MyMPI_Barrier();
        if (rankid == 0)
            cout << "really finished waiting" << endl;
        k = 0;
        ParallelFor
                (Range(ma->GetNFacets()), [&](size_t i) {
                    LocalHeap slh = lh.Split(), &lh = slh;


                    const FacetData &fai = facetdata[i];


                    if (fai.elnr[1] == -1) {
                        auto dp = ma->GetDistantProcs(NodeId(StdNodeType(NT_FACET, ma->GetDimension()), i));
                        if (dp.Size() == 1) {
                            // boundary facet
                            const DGFiniteElement<D> &fel1 = dynamic_cast<const DGFiniteElement<D> &> (fes->GetFE(ElementId(VOL, fai.elnr[0]), lh));
                            const DGFiniteElement<D - 1> &felfacet = dynamic_cast<const DGFiniteElement<D - 1> & > (fes->GetFacetFE(i, lh));

                            IntRange dn1 = fes->GetElementDofs(fai.elnr[0]);

                            size_t ndoffacet = felfacet.GetNDof();
                            size_t ndof1 = fel1.GetNDof();

                            FlatVector<> elu1(ndof1, lh);
                            FlatVector<> trace1(ndoffacet, lh);

                            fel1.GetTrace(fai.facetnr[0], vecu.Range(dn1), trace1);

                            IntegrationRule ir(felfacet.ElementType(), 2 * felfacet.Order());
                            size_t nip = ir.Size();

                            FlatVector<> flown = fai.flown;
                            FlatVector<> tracei1(nip, lh), tracei(nip, lh);

                            felfacet.Evaluate(ir, trace1, tracei1);

                            auto tracei2 = data_get2[dp[0]][k[dp[0]]++];

                            for (size_t j = 0; j < nip; j++)
                                tracei(j) = flown(j) * ((flown(j) > 0) ? tracei1(j) : tracei2[j]);

                            felfacet.EvaluateTrans(ir, tracei, trace1);
                            fel1.GetTraceTrans(fai.facetnr[0], trace1, elu1);

                            {
                                lock_guard<mutex> guard(add_mutex);
                                conv.Range(dn1) += (rankid > dp[0] ? +1 : -1) * elu1;
                            }

                        }
                    }
                });
        MyMPI_Barrier();
    }
};


PYBIND11_MODULE(liblinhyp, m) {

    py::class_<Convection<2>>(m, "Convection")
            .def(py::init<shared_ptr<FESpace>, shared_ptr<CoefficientFunction>>())
            .def("Apply", &Convection<2>::Apply);
}
