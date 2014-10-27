  /*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
  *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
  *
  *  All rights reserved.
  *
  *  This file is part of EMPIRE.
  *
  *  EMPIRE is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  EMPIRE is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
  */
 /***********************************************************************************************//**
 * \file IGAPatchSurfaceTrimming.h
 * This file holds the class IGAPatchSurfaceTrimming.h
 * \date 28/5/2013
 **************************************************************************************************/
 
 #ifndef IGAPatchSurfaceTrimming_H_
 #define IGAPatchSurfaceTrimming_H_
 
 // Inclusion of user defined libraries
 #include "NurbsBasis1D.h"
 #include <vector>
 #include <utility>
 #include <assert.h>

 
 namespace EMPIRE {
     class IGAPatchSurface;
     class IGAControlPoint;
     class Message;
     class IGAPatchSurfaceTrimmingLoop;
     
     /********//**
     * \brief class IGAPatchSurfaceTrimming is a container/proxy to the set of loops trimming the parent Patch
     ***********/
     class IGAPatchSurfaceTrimming {
     private:

     protected:
    	 std::vector<IGAPatchSurfaceTrimmingLoop*> loops;

    	 int outter;

    	 std::vector<std::vector<int> > knotSpanBelonging;
         
         /// The constructor and the destructor and the copy constructor
     public:
         /***********************************************************************************************
          * \brief  Default constructor
          * \author Fabien Pean
          ***********/
         IGAPatchSurfaceTrimming();
         
         /***********************************************************************************************
          * \brief Destructor
          * \author Fabien Pean
          ***********/
         ~IGAPatchSurfaceTrimming();
         /***********************************************************************************************
          * \brief Create data for the knot span matrix
          * \param[in] _uNumKnots The number of knots in U direction
          * \param[in] _vNumKnots The number of knots in V direction
          * \param[in] _knotSpanBelonging The array indicating the knots state, inside,trimmed,outside
          * \author Fabien Pean
          ***********/
         void addTrimInfo(int _uNoKnots, int _vNoKnots, int* _knotSpanBelonging);
         /***********************************************************************************************
          * \brief Setup information about the loop soon to be received
          * \param[in] inner 0 for outter and 1 for inner
          * \param[in] numCurves Number of curves to be received for this loop 
          * \author Fabien Pean
          ***********/
         void addTrimLoop(int _inner, int _numCurves);
         
         /***********************************************************************************************
          * \brief Add a Nurbs curve for the current loop and its attached information
          * \param[in] direction The direction of the curve if is following standard or not
          * \param[in] ID The id of the curve
          * \param[in] _pDegree The polynomial degree of the IGA 1D curve in the u-direction
          * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
          * \param[in] _uKnotVector The underlying knot vector of the IGA 1D curve in the u-direction
          * \param[in] _uNoControlPoints The number of the Control Points for the 2D NURBS patch in the u-direction
          * \param[in] _controlPointNet The set of the Control Points related to the 2D NURBS patch
          * \author Fabien Pean
          ***********/
         void addTrimCurve(int _direction,int _IDBasis, int _pDegree, int _uNoKnots, double* _uKnotVector,
                           int _uNoControlPoints, IGAControlPoint** _controlPointNet);
         
         /***********************************************************************************************
          * \brief Create a linear approximation for every loop using position computed at Greville abscissae
          * 	   And remove non-unique points or aligned points from the set
          * \author Fabien Pean
          ***********/
         void linearizeLoops();
         
         /// Get and set functions
     public:
         /***********************************************************************************************
          * \brief Boolean indicating if there are trimming information
          * \author Fabien Pean
          ***********/
         inline bool isTrimmed() const {
            return (loops.size()>0);
         }
         /***********************************************************************************************
          * \brief Get the outter loop
          * \author Fabien Pean
          ***********/
         inline const IGAPatchSurfaceTrimmingLoop& getOutterLoop() const {
            return *loops[outter];
         }
         /***********************************************************************************************
          * \brief Get the outter loop index
          * \author Fabien Pean
          ***********/
         inline int getOutterLoopIndex() const {
            return outter;
         }
         /***********************************************************************************************
          * \brief Get a specific loop
          * \author Fabien Pean
          ***********/
         inline const IGAPatchSurfaceTrimmingLoop& getLoop(int i) const {
             return *(loops.at(i));
         }
         inline const IGAPatchSurfaceTrimmingLoop& getLastLoop() const {
             return *(loops.back());
         }
         /***********************************************************************************************
          * \brief Get vector of loops
          * \author Fabien Pean
          ***********/
         inline const std::vector<IGAPatchSurfaceTrimmingLoop*>& getLoops() const {
			 return loops;
         }
         /***********************************************************************************************
          * \brief Get the number of loops.
          * \author Fabien Pean
          ***********/
         inline int getNumOfLoops() const {
        	 return loops.size();
         }
         /***********************************************************************************************
		  * \brief Get state of a specific knotSpan
		  * \author Fabien Pean
		  ***********/
		 inline int getKnotSpanInfo(int u, int v) const {
			 return knotSpanBelonging.at(u).at(v);
		 }
		 /***********************************************************************************************
		  * \brief Get the number of loops.
		  * \author Fabien Pean
		  ***********/
		 inline const std::vector<std::vector<int> >& getKnotSpanInfo() const {
			 return knotSpanBelonging;
		 }
     };

     /********//**
     * \brief class IGAPatchSurfaceTrimmingLoop holds all the curves for one trimming
     ***********/
	 class IGAPatchSurfaceTrimmingLoop {
    	 /// The parent class can access it freely
		 friend class IGAPatchSurfaceTrimming;
	 public:
         /***********************************************************************************************
          * \brief Connstructor, reserve data storage for n curves
          * \param[in] _numCurves The number of curves in the loop
          * \author Fabien Pean
          ***********/
		 IGAPatchSurfaceTrimmingLoop(int _numCurves);
         /***********************************************************************************************
          * \brief Destructor
          * \author Fabien Pean
          ***********/
		 ~IGAPatchSurfaceTrimmingLoop(){};
	 private:
         /// The basis functions of the curves
         std::vector<BSplineBasis1D> IGABasis;
         /// Direction = if curve is in the correct orientation C1(1)=C2(0) or not C1(1)=C2(1)
         ///	WARNING : this has not been tested, so code may break if it is used
         std::vector<bool> direction;
         /// Number of control points for each curve
         std::vector<int> uNoControlPoints;
         /// The set of the Control Points of the curves
         std::vector<std::vector<IGAControlPoint*> > ControlPointNet;
         // List of points making up the linearized version of the trimming loop
         std::vector<double> polylines;

         /// Linearizing related functions
     public:
         /***********************************************************************************************
          * \brief Create a linear approximation for this loop
          * \author Fabien Pean
          ***********/
         void linearize();
     private:
         /***********************************************************************************************
          * \brief Compute Greville abscissae for the control point k.
          * The Greville abscissae is the knot where the control point k has the most influence.
          * This should results as the physical shortest distance, or close to it, between curve and the control point
          * \param[in] cp The index of the control point in the local nurbs curve
          * \param[in] pDeg The polynomial degree of the curve
          * \param[in] knotVector The knot vector of the curve
          * \return k The knot vector value corresponding to Greville abscissae
          * \author Fabien Pean
          ***********/
         double computeGrevilleAbscissae(const int cp, const int pDeg, const double*const knotVector);

         void cleanPolygon();

         /// Basis related functions
     public:
         /***********************************************************************************************
          * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters are known
          * \param[in] _curve The index of the curve
          * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the Patch whose surface parameters are _uPrm and _vPrm
          * \param[in] _uPrm The parameter on the u-coordinate line
          * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
          * \author Fabien Pean, Andreas Apostolatos
          ***********/
         void computeCartesianCoordinates(int, double*, double, int);

         /***********************************************************************************************
          * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters and the local basis functions are given
          * \param[in] _curve The index of the curve
          * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm and _vPrm
          * \param[in] _localCoordinates the local coordinates of the point we want to find
          * \compute the knot span Index inside the function. Convenient but in-efficient.
          * \author Fabien Pean, Chenshen Wu
          ***********/
         void computeCartesianCoordinates(int, double*, double);

         /// get functions
	 public:
         /***********************************************************************************************
          * \brief Get the underlying IsoGeometric basis of index i of the loop
          * \author Fabien Pean
          ***********/
         inline const BSplineBasis1D& getIGABasis(int i) const {
             return (const BSplineBasis1D&)IGABasis.at(i);
         }
         /***********************************************************************************************
          * \brief Get the underlying IsoGeometric basis of the loop
          * \author Fabien Pean
          ***********/
         inline const std::vector<BSplineBasis1D>& getIGABasis() const {
             return IGABasis;
         }
         /***********************************************************************************************
          * \brief 	Get the direction for basis i
          * \author Fabien Pean
          ***********/
         inline bool getDirection(int i) const {
             return direction.at(i);
         }
         /***********************************************************************************************
          * \brief Get the number of the Control Points of the curve i
          * \author Fabien Pean
          ***********/
         inline int getNoControlPoints(int i) const {
             return uNoControlPoints.at(i);
         }
         /***********************************************************************************************
          * \brief Get the Control Points of the curve i
          * \author Andreas Apostolatos
          ***********/
         inline const std::vector<IGAControlPoint*>& getControlPointNet(int i) const {
             return ControlPointNet.at(i);
         }
         /***********************************************************************************************
          * \brief Get the linearized version of the curve
          * \author Fabien Pean
          ***********/
         inline const std::vector<double>& getPolylines() const {
             return polylines;
         }
         inline const double* getPolylines(int* size) const {
        	 if(size!=NULL)
        		 *size=polylines.size();
        	 return &polylines[0];
         }
	 };
     
     /***********************************************************************************************
      * \brief Allows for nice debug output
      * \author Fabien Pean
      ***********/
     Message &operator<<(Message &message, const IGAPatchSurfaceTrimming &trim);
     Message &operator<<(Message &message, const IGAPatchSurfaceTrimmingLoop &trim);

     
 }/* namespace EMPIRE */
 
 #endif /* IGAPatchSurfaceTrimming_H_ */
 
