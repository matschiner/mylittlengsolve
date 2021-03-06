#ifndef FILE_MYELEMENT_HPP
#define FILE_MYELEMENT_HPP

/*********************************************************************/
/* File:   myElement.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*
  
My own simple first and second order triangular finite elements

*/

namespace ngfem {

    /*
      A linear triangular finite element is a scalar element in two dimension,
      thus we derive it from the ScalarFiniteElement<2> base class:
     */

    class MyLinearTrig : public ScalarFiniteElement<2> {
    public:
        // constructor
        MyLinearTrig();

        virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

        /*
          Calculate the vector of shape functions in the point ip.
          ip is given in the reference element.
         */
        virtual void CalcShape(const IntegrationPoint &ip,
                               BareSliceVector<> shape) const;

        /*
          Calculate the matrix of derivatives of the shape functions in the point ip.
          dshape is an 3 by 2 matrix in our case.
         */
        virtual void CalcDShape(const IntegrationPoint &ip,
                                BareSliceMatrix<> dshape) const;

        // there are some more functions to bring in ...
        using ScalarFiniteElement<2>::CalcShape;
        using ScalarFiniteElement<2>::CalcDShape;
    };

    class MyLinearRect : public ScalarFiniteElement<2> {
    public:
        // constructor
        MyLinearRect();

        virtual ELEMENT_TYPE ElementType() const { return ET_QUAD; }

        /*
          Calculate the vector of shape functions in the point ip.
          ip is given in the reference element.
         */
        virtual void CalcShape(const IntegrationPoint &ip,
                               BareSliceVector<> shape) const;

        /*
          Calculate the matrix of derivatives of the shape functions in the point ip.
          dshape is an 3 by 2 matrix in our case.
         */
        virtual void CalcDShape(const IntegrationPoint &ip,
                                BareSliceMatrix<> dshape) const;

        // there are some more functions to bring in ...
        using ScalarFiniteElement<2>::CalcShape;
        using ScalarFiniteElement<2>::CalcDShape;
    };

    class MyCubicTrig : public ScalarFiniteElement<2> {
    public:
        // constructor
        MyCubicTrig();

        virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }


        virtual void CalcShape(const IntegrationPoint &ip,
                               BareSliceVector<> shape) const;

        virtual void CalcDShape(const IntegrationPoint &ip,
                                BareSliceMatrix<> dshape) const;

        virtual std::array<AutoDiff<2>, 10> GetBasis(const IntegrationPoint &ip) const;

        // there are some more functions to bring in ...
        using ScalarFiniteElement<2>::CalcShape;
        using ScalarFiniteElement<2>::CalcDShape;
    };
    class ThirdOrderLineSegment : public ScalarFiniteElement<1>
    {
    public:
        ThirdOrderLineSegment();
        virtual ELEMENT_TYPE ElementType() const { return ET_SEGM; }
        virtual void CalcShape(const IntegrationPoint &ip, BareSliceVector<> shape) const;
        virtual void CalcDShape(const IntegrationPoint &ip, BareSliceMatrix<> dshape) const;

        using ScalarFiniteElement<1>::CalcShape;
        using ScalarFiniteElement<1>::CalcDShape;
    };
    /*
      A triangular finite element with second order basis functions
     */
    class MyQuadraticTrig : public ScalarFiniteElement<2> {
    public:
        MyQuadraticTrig();

        virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

        virtual void CalcShape(const IntegrationPoint &ip,
                               BareSliceVector<> shape) const;

        virtual void CalcDShape(const IntegrationPoint &ip,
                                BareSliceMatrix<> dshape) const;

        using ScalarFiniteElement<2>::CalcShape;
        using ScalarFiniteElement<2>::CalcDShape;
    };

}

void ExportMyElement(py::module m);

#endif // FILE_MYELEMENT_HPP

