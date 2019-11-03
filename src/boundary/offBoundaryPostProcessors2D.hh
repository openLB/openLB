/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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

#ifndef OFF_BOUNDARY_POST_PROCESSORS_2D_HH
#define OFF_BOUNDARY_POST_PROCESSORS_2D_HH

#include "offBoundaryPostProcessors2D.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "core/cell.h"

namespace olb {

/////////// LinearBouzidiPostProcessor2D /////////////////////////////////////

/* Bouzidi Interpolation scheme of first order
 *
 * fluid nodes               wall  solid node
 * --o-------<-o->-----<-o->--|----x----
 *            xB         x        xN
 * directions: --> iPop
 *             <-- opp
 *
*/

template<typename T, typename DESCRIPTOR>
ZeroVelocityBouzidiLinearPostProcessor2D<T,DESCRIPTOR>::
ZeroVelocityBouzidiLinearPostProcessor2D(int x_, int y_, int iPop_, T dist_)
  : x(x_), y(y_), iPop(iPop_), dist(dist_)
{
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  opp = util::opposite<L>(iPop);
  xN = x + descriptors::c<L>(iPop,0);
  yN = y + descriptors::c<L>(iPop,1);

  if (dist >= 0.5) {
    xB = x - descriptors::c<L>(iPop,0);
    yB = y - descriptors::c<L>(iPop,1);
    q = 1/(2*dist);
    iPop2 = opp;
  } else {
    xB = x;
    yB = y;
    q = 2*dist;
    iPop2 = iPop;
  }
  /*
    std::cout << "ZeroVelocityLinear (" << x << "," << y << "," <<
      "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
      "), opp: " << opp << ", bP: (" << xB << "," << yB << "," <<
      "), dist: " << dist << ", q: " << q << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBouzidiLinearPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBouzidiLinearPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  blockLattice.get(x, y)[opp] = q*blockLattice.get(xN, yN)[iPop] +
                                (1-q)*blockLattice.get(xB, yB)[iPop2];
}

template<typename T, typename DESCRIPTOR>
VelocityBouzidiLinearPostProcessor2D<T,DESCRIPTOR>::
VelocityBouzidiLinearPostProcessor2D(int x_, int y_, int iPop_, T dist_)
  : x(x_), y(y_), iPop(iPop_), dist(dist_)
{
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  opp = util::opposite<L>(iPop);
  xN = x + descriptors::c<L>(iPop,0);
  yN = y + descriptors::c<L>(iPop,1);

  if (dist >= 0.5) {
    xB = x - descriptors::c<L>(iPop,0);
    yB = y - descriptors::c<L>(iPop,1);
    q = 1/(2*dist);
    ufrac = q;
    iPop2 = opp;
  } else {
    xB = x;
    yB = y;
    q = 2*dist;
    iPop2 = iPop;
    ufrac = 1;
  }
  /*
    std::cout << "VelocityLinear (" << x << "," << y << "," <<
      "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
      "), opp: " << opp << ", bP: (" << xB << "," << yB << "," <<
      "), dist: " << dist << ", q: " << q << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void VelocityBouzidiLinearPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void VelocityBouzidiLinearPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(xN, yN);
  T velCoeff = ufrac*dynamics->getVelocityCoefficient(iPop);
  dynamics->defineRho( blockLattice.get(xN, yN), blockLattice.get(x, y).computeRho() );
  blockLattice.get(x, y)[opp] = q*blockLattice.get(xN, yN)[iPop] +
                                (1-q)*blockLattice.get(xB, yB)[iPop2] + velCoeff;
}


//////// CornerBouzidiPostProcessor2D ///////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBounceBackPostProcessor2D<T,DESCRIPTOR>::
ZeroVelocityBounceBackPostProcessor2D(int x_, int y_, int iPop_, T dist_)
  : x(x_), y(y_), iPop(iPop_), dist(dist_)
{
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  opp = util::opposite<L>(iPop);
  xN = x + descriptors::c<L>(iPop,0);
  yN = y + descriptors::c<L>(iPop,1);
  /*
    std::cout << "Corner (" << x << "," << y << "," <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBounceBackPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBounceBackPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  blockLattice.get(x, y)[opp] = blockLattice.get(xN, yN)[iPop];
}


template<typename T, typename DESCRIPTOR>
VelocityBounceBackPostProcessor2D<T,DESCRIPTOR>::
VelocityBounceBackPostProcessor2D(int x_, int y_, int iPop_, T dist_)
  : x(x_), y(y_), iPop(iPop_), dist(dist_)
{
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  opp = util::opposite<L>(iPop);
  xN = x + descriptors::c<L>(iPop,0);
  yN = y + descriptors::c<L>(iPop,1);

  /*
    std::cout << "Corner (" << x << "," << y << "," <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void VelocityBounceBackPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void VelocityBounceBackPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(xN, yN);
  T velCoeff = dynamics->getVelocityCoefficient(iPop);
  dynamics->defineRho( blockLattice.get(xN, yN), blockLattice.get(x, y).computeRho() );
  blockLattice.get(x, y)[opp] = blockLattice.get(xN, yN)[iPop] + velCoeff;
}


template<typename T, typename DESCRIPTOR>
AntiBounceBackPostProcessor2D<T,DESCRIPTOR>::
AntiBounceBackPostProcessor2D(int x_, int y_, int iPop_)
  : x(x_), y(y_), iPop(iPop_)
{
  typedef DESCRIPTOR L;
  opp = util::opposite<L>(iPop);
  xN = x + descriptors::c<L>(iPop,0);
  yN = y + descriptors::c<L>(iPop,1);

  /*
    std::cout << "Corner (" << x << "," << y << "," <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void AntiBounceBackPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void AntiBounceBackPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  /*Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(xN, yN);
  T velCoeff = dynamics->getVelocityCoefficient(iPop);
  dynamics->defineRho( blockLattice.get(xN, yN), blockLattice.get(x, y).computeRho() );*/
  if (descriptors::c<DESCRIPTOR>(iPop,1)==0) {
    blockLattice.get(x, y)[opp] = -blockLattice.get(xN, yN)[iPop];  // + velCoeff;
  }
  //std::cout << "here" << std::endl;
}


template<typename T, typename DESCRIPTOR>
BoundaryStreamPostProcessor2D<T,DESCRIPTOR>::
BoundaryStreamPostProcessor2D(int x_, int y_, const bool streamDirection[DESCRIPTOR::q])
  : x(x_), y(y_)
{
  for (int iPop = 0; iPop < DESCRIPTOR::q ; ++iPop) {
    this->_streamDirections[iPop] = streamDirection[iPop];
  }
}

template<typename T, typename DESCRIPTOR>
void BoundaryStreamPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void BoundaryStreamPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    if (_streamDirections[iPop]) {
      blockLattice.get(x + descriptors::c<DESCRIPTOR>(iPop,0), y + descriptors::c<DESCRIPTOR>(iPop,1))[iPop] = blockLattice.get(x, y)[iPop];
    }
  }
}


////////  LinearBouzidiBoundaryPostProcessorGenerator ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBouzidiLinearPostProcessorGenerator2D<T,DESCRIPTOR>::
ZeroVelocityBouzidiLinearPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
ZeroVelocityBouzidiLinearPostProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new ZeroVelocityBouzidiLinearPostProcessor2D<T,DESCRIPTOR>
         ( this->x, this->y, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
ZeroVelocityBouzidiLinearPostProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator2D<T,DESCRIPTOR>
         (this->x, this->y, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
VelocityBouzidiLinearPostProcessorGenerator2D<T,DESCRIPTOR>::
VelocityBouzidiLinearPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
VelocityBouzidiLinearPostProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new VelocityBouzidiLinearPostProcessor2D<T,DESCRIPTOR>
         ( this->x, this->y, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
VelocityBouzidiLinearPostProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new VelocityBouzidiLinearPostProcessorGenerator2D<T,DESCRIPTOR>
         (this->x, this->y, this->iPop, this->dist);
}

/////////// CornerBouzidiBoundaryPostProcessorGenerator /////////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>::
ZeroVelocityBounceBackPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
ZeroVelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new ZeroVelocityBounceBackPostProcessor2D<T,DESCRIPTOR>
         ( this->x, this->y, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
ZeroVelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new ZeroVelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>
         (this->x, this->y, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
VelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>::
VelocityBounceBackPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
VelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new VelocityBounceBackPostProcessor2D<T,DESCRIPTOR>
         ( this->x, this->y, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
VelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new VelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>
         (this->x, this->y, this->iPop, this->dist);
}


template<typename T, typename DESCRIPTOR>
AntiBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>::
AntiBounceBackPostProcessorGenerator2D(int x_, int y_, int iPop_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
AntiBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new AntiBounceBackPostProcessor2D<T,DESCRIPTOR>
         ( this->x, this->y, this->iPop);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
AntiBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new AntiBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>
         (this->x, this->y, this->iPop);
}

template<typename T, typename DESCRIPTOR>
BoundaryStreamPostProcessorGenerator2D<T,DESCRIPTOR>::
BoundaryStreamPostProcessorGenerator2D(int x_, int y_, const bool streamDirections[DESCRIPTOR::q])
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x_, x_, y_, y_),
    x(x_), y(y_)
{
  for (int iPop = 0; iPop < DESCRIPTOR::q ; ++iPop) {
    this->_streamDirections[iPop] = streamDirections[iPop];
  }
}

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
BoundaryStreamPostProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new BoundaryStreamPostProcessor2D<T,DESCRIPTOR>
         ( this->x, this->y, this->_streamDirections);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
BoundaryStreamPostProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new BoundaryStreamPostProcessorGenerator2D<T,DESCRIPTOR>
         (this->x, this->y, this->_streamDirections);
}
}  // namespace olb

#endif
