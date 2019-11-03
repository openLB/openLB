/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, Orestis Malaspinas and Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

/** \file
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef INAMURO_BOUNDARY_2D_HH
#define INAMURO_BOUNDARY_2D_HH

#include "inamuroBoundary2D.h"
#include "inamuroAnalyticalDynamics.h"
#include "inamuroAnalyticalDynamics.hh"
#include "boundaryInstantiator2D.h"

namespace olb {

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class InamuroBoundaryManager2D {
public:
  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv=NULL);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getExternalVelocityCornerProcessor(int x, int y);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getInternalVelocityCornerProcessor(int x, int y);
};

////////// InamuroBoundaryManager2D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR, VelocityBM, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new InamuroAnalyticalDynamics<T,DESCRIPTOR, MixinDynamics, direction, orientation>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR, PressureBM, direction,orientation>;
}
template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new InamuroAnalyticalDynamics<T,DESCRIPTOR, MixinDynamics, direction, orientation>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y)
{
  return new OuterVelocityCornerProcessorGenerator2D<T,DESCRIPTOR, xNormal,yNormal> (x,y);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM2D<T,DESCRIPTOR, xNormal,yNormal>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
InamuroBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerProcessor
(int x, int y)
{
  return nullptr;
}

////////// Factory function //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeBoundaryCondition2D<T,DESCRIPTOR>* createInamuroBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block)
{
  return new BoundaryConditionInstantiator2D <
         T, DESCRIPTOR,
         InamuroBoundaryManager2D<T,DESCRIPTOR, MixinDynamics> > (block);
}

}  // namespace olb

#endif
