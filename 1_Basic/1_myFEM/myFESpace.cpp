/*********************************************************************/
/* File:   myFESpace.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/


/*

My own FESpace for linear and quadratic triangular elements.

A fe-space provides the connection between the local reference
element, and the global mesh.

*/


#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>
#include <python_ngstd.hpp>
#include "myElement.hpp"
#include "myFESpace.hpp"


namespace ngcomp {

    MyFESpace::MyFESpace(shared_ptr<MeshAccess> ama, const Flags &flags)
            : FESpace(ama, flags) {
        cout << "Constructor of MyFESpace" << endl;
        //cout << "Flags = " << flags << endl;

        order = flags.GetNumFlag("order", 1);
        FE_geom = flags.GetStringFlag("FE_geom", "trig");
        cout << "You have chosen the following type of FE elements: " << FE_geom << endl;

        cout << "You have chosen the following order of elements: " << order << endl;

        // needed for symbolic integrators and to draw solution
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();

        // (still) needed to draw solution
        integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(),
                                                     make_shared<ConstantCoefficientFunction>(1));
    }

    DocInfo MyFESpace::GetDocu() {
        auto docu = FESpace::GetDocu();
        docu.Arg("secondorder") = "bool = False\n"
                                  "  Use second order basis functions";
        docu.Arg("order") = "order = 1,2,3 \n"
                            "  Use higher order basis functions";
        return docu;
    }

    void MyFESpace::Update(LocalHeap &lh) {
        // some global update:
        cout << "Update MyFESpace, #vert = " << ma->GetNV()
             << ", #edge = " << ma->GetNEdges() << ", #elems = " << ma->GetNElements(2) << endl;

        // number of vertices
        nvert = ma->GetNV();
        nedges = ma->GetNEdges();
        nelems = ma->GetNE();

        // number of dofs:
        ndof = nvert;
        if (order == 2) {
            ndof += nedges;  // num vertics + num edges
        } else if (order == 3) {
            ndof += 2 * nedges;  // num vertics + num edges
            ndof += nelems;

        }
    }

    void MyFESpace::GetDofNrs(ElementId ei, Array<DofId> &dnums) const {
        // returns dofs of element ei
        // may be a volume triangle or boundary segment


        dnums.SetSize(0);

        // first dofs are vertex numbers:
        if (order <= 2) {
            for (auto v : ma->GetElVertices(ei))
                dnums.Append(v);
            if (order == 2) {
                // more dofs on edges:
                for (auto e : ma->GetElEdges(ei))
                    dnums.Append(nvert + e);
            }
        }


        if (order == 3) {
            for (auto v : ma->GetElVertices(ei)) {
                dnums.Append(v);
            }

            Array<int> nb;

            for (auto e :  ma->GetElEdges(ei)) {
                ma->GetEdgeElements(e, nb);
                if (ei.Nr() <= nb[0] && (nb.Size() == 1 || ei.Nr() <= nb[1])) {
                    dnums.Append(nvert + 2 * e);
                    dnums.Append(nvert + 2 * e + 1);
                } else {
                    dnums.Append(nvert + 2 * e + 1);
                    dnums.Append(nvert + 2 * e);
                }
            }
            if (ei.IsVolume()) {
                dnums.Append(nvert + 2 * nedges + ei.Nr());
            }

        }

    }

    FiniteElement &MyFESpace::GetFE(ElementId ei, Allocator &alloc) const {
        if (ei.IsVolume()) {


            if (order == 1) {
                if (ma->GetElVertices(ei).Size() == 4) {
                    return *new(alloc) MyLinearRect;
                } else if (ma->GetElVertices(ei).Size() == 3) {
                    return *new(alloc) MyLinearTrig;
                } else {
                    cout << "no method found for FE with " << ma->GetElVertices(ei).Size() << " vertices" << endl;
                    return *new(alloc) MyLinearTrig;
                }
            } else if (order == 2) {
                return *new(alloc) MyQuadraticTrig;
            } else if (order == 3) {
                return *new(alloc) MyCubicTrig;
            }

        } else {
            if (order == 1)
                return *new(alloc) FE_Segm1;
            else if (order == 3)
                return *new(alloc) ThirdOrderLineSegment;
            else {
                cout << "not line Segment specified" << endl;
                return *new(alloc) FE_Segm2;
            }
        }
    }

    /*
      register fe-spaces
      Object of type MyFESpace can be defined in the pde-file via
      "define fespace v -type=myfespace"
    */

    static RegisterFESpace<MyFESpace> initifes("myfespace");
}

void ExportMyFESpace(py::module m) {
    using namespace ngcomp;
    /*
      We just export the class here and use the FESpace constructor to create our space.
      This has the advantage, that we do not need to specify all the flags to parse (like
      dirichlet, definedon,...), but we can still append new functions only for that space.
    */
    auto docu = MyFESpace::GetDocu();
    auto myfes = py::class_<MyFESpace, shared_ptr<MyFESpace>, FESpace>
            (m, "MyFESpace", (docu.short_docu + "\n\n" + docu.long_docu).c_str());
    myfes
            /*
               this is optional, if you don't write an init function, you can create your fespace
               with FESpace("myfes",mesh,...), but it's nicer to write MyFESpace(mesh,...) ;)
            */
            .def(py::init([myfes](shared_ptr<MeshAccess> ma, py::kwargs kwa) {
                py::list info;
                info.append(ma);
                auto flags = CreateFlagsFromKwArgs(myfes, kwa, info);
                auto fes = make_shared<MyFESpace>(ma, flags);
                LocalHeap glh(100000000, "init-fes-lh");
                fes->Update(glh);
                fes->FinalizeUpdate(glh);
                return fes;
            }), py::arg("mesh"))
                    /*
                      this is, so that we do not get an 'undocumented flag' warning
                    */
            .def_static("__flags_doc__", [docu]() {
                auto doc = py::cast<py::dict>(py::module::import("ngsolve").
                        attr("FESpace").
                        attr("__flags_doc__")());
                for (auto &flagdoc : docu.arguments)
                    doc[get<0>(flagdoc).c_str()] = get<1>(flagdoc);
                return doc;
            })
            .def("GetNVert", &MyFESpace::GetNVert);
}
