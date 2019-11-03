/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Drag force models for Lagrangian two-way coupling -- generic implementation.
 */

#ifndef LB_DRAG_MODELS_HH
#define LB_DRAG_MODELS_HH

namespace olb {


////////////////////// Class DragModelBase ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
DragModelBase<T,Lattice,Particle>::DragModelBase(UnitConverter<T, Lattice>& converter)
         : _converter(converter)
{}


////////////////////// Class StokesSimplifiedDragModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
StokesSimplifiedDragModel<T,Lattice,Particle>::StokesSimplifiedDragModel(UnitConverter<T, Lattice>& converter)
         : DragModelBase<T,Lattice,Particle>(converter)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
T StokesSimplifiedDragModel<T,Lattice,Particle>::operator() (
                         Particle<T>* p, T latticeVelF[], T latticeVelP[], int globicFull[] )
{
  return 1.83;
}


////////////////////// Class MorsiDragModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
MorsiDragModel<T,Lattice,Particle>::MorsiDragModel(UnitConverter<T, Lattice>& converter)
         : DragModelBase<T,Lattice,Particle>(converter)
{
  this->_reP = std::make_shared<NewtonianParticleReynoldsNumber<T,Lattice,Particle> > (this->_converter);
}

template<typename T, typename Lattice, template<typename V> class Particle>
T MorsiDragModel<T,Lattice,Particle>::operator() (
                         Particle<T>* p, T latticeVelF[], T latticeVelP[], int globicFull[] )
{
  T physVelRelative = this->_converter.getPhysVelocity (
                             sqrt( pow(latticeVelF[0] - latticeVelP[0],2) +
                                   pow(latticeVelF[1] - latticeVelP[1],2) +
                                   pow(latticeVelF[2] - latticeVelP[2],2) ) );

  T ReP = this->_reP->operator() (p, physVelRelative, globicFull);

  T a[3] = {T(), T(), T()};
  if (ReP < 0.1) {
    a[0] = 0.0;     a[1] = 24.0;      a[2] = 0.0;
  }
  else if (ReP < 1.0) {
    a[0] = 3.69;    a[1] = 22.73;     a[2] = 0.0903;
  }
  else if (ReP < 10.0) {
    a[0] = 1.222;   a[1] = 29.16667;  a[2] =-3.8889;
  }
  else if (ReP < 100.0) {
    a[0] = 0.6167;  a[1] = 46.5;      a[2] =-116.67;
  }
  else if (ReP < 1000.0) {
    a[0] = 0.3644;  a[1] = 498.33;    a[2] =-2778;
  }
  else if (ReP < 5000.0) {
    a[0] = 0.357;   a[1] = 148.62;    a[2] =-4.75e4;
  }
  else if (ReP < 10000.0) {
    a[0] = 0.46;    a[1] =-490.546;   a[2] = 57.87e4;
  }
  else {
    a[0] = 0.5191;  a[1] =-1662.5;    a[2] = 5.4167e6;
  }


  return ( a[0] + a[1]/ReP + a[2]/(ReP*ReP) ) * physVelRelative;
}


////////////////////// Class DewsburyDragModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
DewsburyDragModel<T,Lattice,Particle>::DewsburyDragModel(UnitConverter<T, Lattice>& converter)
         : DragModelBase<T,Lattice,Particle>(converter)
{
  this->_reP = std::make_shared<NewtonianParticleReynoldsNumber<T,Lattice,Particle> > (this->_converter);
}

template<typename T, typename Lattice, template<typename V> class Particle>
T DewsburyDragModel<T,Lattice,Particle>::operator() (
                         Particle<T>* p, T latticeVelF[], T latticeVelP[], int globicFull[] )
{
  T physVelRelative = this->_converter.getPhysVelocity (
                             sqrt( pow(latticeVelF[0] - latticeVelP[0],2) +
                                   pow(latticeVelF[1] - latticeVelP[1],2) +
                                   pow(latticeVelF[2] - latticeVelP[2],2) ) );

  T ReP = this->_reP->operator() (p, physVelRelative, globicFull);

  T Cd = 0.95;
  if (ReP <= 195.) {
    Cd = 16./ReP * (1. + 0.173*pow(ReP, 0.657))
       + 0.413 / (1. + 16300*pow(ReP, -1.09));
  }
  return Cd * this->_converter.getLatticeVelocity(physVelRelative);
}


////////////////////// Class PowerLawDewsburyDragModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
PowerLawDewsburyDragModel<T,Lattice,Particle>::PowerLawDewsburyDragModel (
        UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice )
         : DewsburyDragModel<T,Lattice,Particle>(converter)
{
  this->_reP = std::make_shared<PowerLawParticleReynoldsNumber<T,Lattice,Particle> > (this->_converter, sLattice);
}


}

#endif
