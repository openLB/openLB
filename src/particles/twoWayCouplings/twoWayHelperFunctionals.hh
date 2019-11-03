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

/* Helper functionals for Lagrangian two-way coupling methods -- generic implementation.
 */

#ifndef LB_TWO_WAY_HELPER_FUNCTIONALS_HH
#define LB_TWO_WAY_HELPER_FUNCTIONALS_HH

namespace olb {

////////////////////// Class ParticleReynoldsNumberBase ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
ParticleReynoldsNumberBase<T,Lattice,Particle>::ParticleReynoldsNumberBase(UnitConverter<T, Lattice>& converter)
         : _converter(converter)
{}


////////////////////// Class NewtonianParticleReynoldsNumber ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NewtonianParticleReynoldsNumber<T,Lattice,Particle>::NewtonianParticleReynoldsNumber(UnitConverter<T, Lattice>& converter)
         : ParticleReynoldsNumberBase<T,Lattice,Particle>(converter)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
T NewtonianParticleReynoldsNumber<T,Lattice,Particle>::operator() ( Particle<T>* p, T magU, int globicFull[])
{
  T ReP = 2. * p->getRad() * magU / this->_converter.getPhysViscosity();
  return ReP > this->_RePmin ? ReP : this->_RePmin;
}


////////////////////// Class PowerLawParticleReynoldsNumber ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
PowerLawParticleReynoldsNumber<T,Lattice,Particle>::PowerLawParticleReynoldsNumber (
        UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice )
         : ParticleReynoldsNumberBase<T,Lattice,Particle>(converter),
           _sLattice(sLattice)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
T PowerLawParticleReynoldsNumber<T,Lattice,Particle>::operator() ( Particle<T>* p, T magU, int globicFull[])
{
  // loc() indicates the local cuboid number locIC of the actual processing thread,
  // for given global cuboid number iC
  // this is to get appropriate particle system on locIC
  int locIC = _sLattice.getLoadBalancer().loc(globicFull[0]);

  T physPosP[3] = {T(), T(), T()}; // particle's physical position
  physPosP[0] = (p->getPos()[0]);
  physPosP[1] = (p->getPos()[1]);
  physPosP[2] = (p->getPos()[2]);

  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = { globicFull[1], globicFull[2], globicFull[3] };

  // getting the power-law relaxation frequency form the dynamics's external field
  T omega = _sLattice.getBlockLattice(locIC).get ( 
                      latticeRoundedPosP[0],
                      latticeRoundedPosP[1],
                      latticeRoundedPosP[2] ).template getFieldPointer<descriptors::OMEGA>()[0];

  // physical viscosity from relaxation time
  T nu = this->_converter.getPhysViscosity (
                          (1./omega - 0.5) / descriptors::invCs2<T,Lattice>() );
  
  T ReP = 2. * p->getRad() * magU / nu;
  return ReP > this->_RePmin ? ReP : this->_RePmin;
}


////////////////////// Class TwoWayHelperFunctional ////////////////////////

template<typename T, typename Lattice>
TwoWayHelperFunctional<T, Lattice>::TwoWayHelperFunctional (
                     UnitConverter<T, Lattice>& converter,
                     SuperLattice3D<T, Lattice>& sLattice )
         : _converter(converter),
           _sLattice(sLattice)
{}

template<typename T, typename Lattice>
TwoWayHelperFunctional<T, Lattice>::~TwoWayHelperFunctional()
{
  _interpLatticeDensity.reset();
  _interpLatticeVelocity.reset();
}


////////////////////// Class NaiveMomentumExchange ////////////////////////

template<typename T, typename Lattice>
NaiveMomentumExchange<T, Lattice>::NaiveMomentumExchange (
                     UnitConverter<T, Lattice>& converter,
                     SuperLattice3D<T, Lattice>& sLattice,
                     std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > interpLatticeDensity )
         : TwoWayHelperFunctional<T, Lattice>(converter, sLattice)
{
  this->_interpLatticeDensity = interpLatticeDensity;
}

template<typename T, typename Lattice>
bool NaiveMomentumExchange<T, Lattice>::operator() ( T gF[], T latticeVelF[], T latticeVelP[],
                            T physPosP[], int latticeRoundedP[],
                            int globic )
{
  T magU = sqrt( pow(latticeVelF[0] - latticeVelP[0],2) +
                 pow(latticeVelF[1] - latticeVelP[1],2) +
                 pow(latticeVelF[2] - latticeVelP[2],2) );

  // Interpolated/exact pdf (depending on the functional)
  T f[Lattice::q] = { T() };
  this->_interpLatticeDensity->operator()(f, physPosP, globic);

  // Mock cell with the same dynamics of the cell containing the particle
  int locIC = this->_sLattice.getLoadBalancer().loc(globic);
  Cell<T,Lattice> cell ( this->_sLattice.getBlockLattice(locIC).get ( 
                         latticeRoundedP[0],
                         latticeRoundedP[1],
                         latticeRoundedP[2] ).getDynamics() );

  // Filling the mock cell with the interpolated pdf
  for (unsigned iPop = 0; iPop < Lattice::q; ++iPop) {
    cell[iPop] = f[iPop];
  }

  T rhoL = cell.computeRho();

  gF[0] = rhoL * magU;
  gF[1] = rhoL * magU;
  gF[2] = rhoL * magU;

  return true;
}


////////////////////// Class LaddMomentumExchange ////////////////////////

template<typename T, typename Lattice>
LaddMomentumExchange<T, Lattice>::LaddMomentumExchange (
                     UnitConverter<T, Lattice>& converter,
                     SuperLattice3D<T, Lattice>& sLattice,
                     std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > interpLatticeDensity,
                     std::shared_ptr<SuperLatticeInterpPhysVelocity3D<T, Lattice> > interpLatticeVelocity )
         : TwoWayHelperFunctional<T, Lattice>(converter, sLattice)
{
  this->_interpLatticeDensity = interpLatticeDensity;
  this->_interpLatticeVelocity = interpLatticeVelocity;
}

template<typename T, typename Lattice>
bool LaddMomentumExchange<T, Lattice>::operator() ( T gF[], T latticeVelF[], T latticeVelP[],
                            T physPosP[], int latticeRoundedP[],
                            int globic )
{
  T physLatticeL = this->_converter.getConversionFactorLength();

  // force density gF
  gF[0] = T();
  gF[1] = T();
  gF[2] = T();

  T fiPop = T();
  T sp = T(); // dot product for particle velocity
  T faPos[3] = {T(), T(), T()}; // fAlphaPosition = particle position
  T fbPos[3] = {T(), T(), T()}; // fBetaPosition = neighbor position to particle position in direction iPop

  T fa[Lattice::q] = { T() }; // fAlpha = interpolated density to fAlphaPosition
  T fb[Lattice::q] = { T() }; // fBeta = interpolated density to fBetaPosition
  T lFU[3] = {T(), T(), T()};

  // runs through all q discrete velocity directions
  for (unsigned iPop = 0; iPop < Lattice::q; ++iPop) {
    // physical position on neighbor to get pre-streaming collision part
    faPos[0] = physPosP[0] + physLatticeL * descriptors::c<Lattice>(iPop,0);
    faPos[1] = physPosP[1] + physLatticeL * descriptors::c<Lattice>(iPop,1);
    faPos[2] = physPosP[2] + physLatticeL * descriptors::c<Lattice>(iPop,2);
    // Lagrange interpolated polynomial to get density on particle position
    this->_interpLatticeDensity->operator() (fa, faPos, globic);

    // physical position on neighbor to get pre-streaming collision part
    fbPos[0] = physPosP[0] - physLatticeL * descriptors::c<Lattice>(iPop,0);
    fbPos[1] = physPosP[1] - physLatticeL * descriptors::c<Lattice>(iPop,1);
    fbPos[2] = physPosP[2] - physLatticeL * descriptors::c<Lattice>(iPop,2);
    // Lagrange interpolated polynomial to get density on particle position
    this->_interpLatticeDensity->operator() (fb, fbPos, globic);

    // fiPop = density on fBetaPosition in direction iPop
    fiPop = fb[util::opposite<Lattice >(iPop)];
    // Get f_l of the boundary cell
    // add density fAlphaL of opposite direction to iPop
    fiPop -= fa[iPop];

    // physical velocity
    lFU[0] = -descriptors::c<Lattice>(iPop,0) * fiPop;
    lFU[1] = -descriptors::c<Lattice>(iPop,1) * fiPop;
    lFU[2] = -descriptors::c<Lattice>(iPop,2) * fiPop;

    // point product
    sp = descriptors::c<Lattice>(iPop,0) * latticeVelP[0] + descriptors::c<Lattice>(iPop,1) * latticeVelP[1]
        + descriptors::c<Lattice>(iPop,2) * latticeVelP[2];
    sp *= 2. * descriptors::invCs2<T,Lattice>()
          * descriptors::t<T,Lattice>(iPop);

    // external force density that acts on particles
    gF[0] += (lFU[0] - descriptors::c<Lattice>(iPop,0) * (sp));
    gF[1] += (lFU[1] - descriptors::c<Lattice>(iPop,1) * (sp));
    gF[2] += (lFU[2] - descriptors::c<Lattice>(iPop,2) * (sp));
  }
  gF[0] = fabs(gF[0]);
  gF[1] = fabs(gF[1]);
  gF[2] = fabs(gF[2]);

  return true;
}


////////////////////// Class SmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
SmoothingFunctional<T, Lattice>::SmoothingFunctional (
                 T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice )
         : _kernelLength(kernelLength),
           _converter(converter),
           _sLattice(sLattice)
{}

template<typename T, typename Lattice>
int SmoothingFunctional<T, Lattice>::getSize()
{
  return _latticePosAndWeight.size();
}

template<typename T, typename Lattice>
void SmoothingFunctional<T, Lattice>::getLatticePos(int latticePos[], int i)
{
  latticePos[0] = _latticePosAndWeight[i].latticePos[0];
  latticePos[1] = _latticePosAndWeight[i].latticePos[1];
  latticePos[2] = _latticePosAndWeight[i].latticePos[2];
}

template<typename T, typename Lattice>
int SmoothingFunctional<T, Lattice>::getGlobic(int i)
{
  return _latticePosAndWeight[i].globic;
}

template<typename T, typename Lattice>
T SmoothingFunctional<T, Lattice>::getWeight(int i)
{
  return _latticePosAndWeight[i].weight;
}

template<typename T, typename Lattice>
bool SmoothingFunctional<T, Lattice>::update(T physPosP[], int globic)
{
  // Bottom-left corner of a cube centered at the particle, with side 2*_kernelLength
  T physPosMin[3] = {T(), T(), T()};
  physPosMin[0] = physPosP[0] - _kernelLength;
  physPosMin[1] = physPosP[1] - _kernelLength;
  physPosMin[2] = physPosP[2] - _kernelLength;
  int latticePosMin[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
           latticePosMin, physPosMin );

  // Top-right corner of a cube centered at the particle, with side 2*_kernelLength
  T physPosMax[3] = {T(), T(), T()};
  physPosMax[0] = physPosP[0] + _kernelLength;
  physPosMax[1] = physPosP[1] + _kernelLength;
  physPosMax[2] = physPosP[2] + _kernelLength;
  int latticePosMax[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
           latticePosMax, physPosMax );

  // Clearing the _latticePosAndWeight list
  _latticePosAndWeight.clear();

  T normalizer = T();
  int iLatticePos[3] = {0, 0, 0};
  // Cycling all the cells on a cube containing a sphee centered in bubble's position and with kernel smoothing length as radius
  for (iLatticePos[0]=latticePosMin[0]; iLatticePos[0]<=latticePosMax[0]; iLatticePos[0]++) {
    for (iLatticePos[1]=latticePosMin[1]; iLatticePos[1]<=latticePosMax[1]; iLatticePos[1]++) {
      for (iLatticePos[2]=latticePosMin[2]; iLatticePos[2]<=latticePosMax[2]; iLatticePos[2]++) {

        T iPhysPos[3] = {T(), T(), T()};
        this->_sLattice.getCuboidGeometry().get(globic).getPhysR (
               iPhysPos, iLatticePos );

        // Is the voxel within a smooting kernel length from the bubble's position?
        if ( pow(physPosP[0] - iPhysPos[0], 2) +
             pow(physPosP[1] - iPhysPos[1], 2) +
             pow(physPosP[2] - iPhysPos[2], 2) < pow(_kernelLength, 2) ) {

          // Adding the voxel's position (and relative weight) to the _latticePosAndWeight list
          LatticePosAndWeight<T> item;
          item.latticePos[0] = iLatticePos[0];
          item.latticePos[1] = iLatticePos[1];
          item.latticePos[2] = iLatticePos[2];
          item.weight = this->compute(physPosP, iPhysPos);

          normalizer += item.weight;
          _latticePosAndWeight.push_back(item);
        }
      }
    }
  }

  // If normalizer is zero, then no voxels are within a kernel smoothing length from the bubble's location.
  // And it is a problem.
  if (normalizer == T()) {
    std::cout << "ERROR: SmoothingFunctional::update(...):" << std::endl
              << "[smoothingFunctional] physPosP: "
              << physPosP[0] << " "
              << physPosP[1] << " "
              << physPosP[2] << std::endl
              << "[smoothingFunctional] physPosMin: "
              << physPosMin[0] << " "
              << physPosMin[1] << " "
              << physPosMin[2] << std::endl
              << "[smoothingFunctional] physPosMax: "
              << latticePosMax[0] << " "
              << latticePosMax[1] << " "
              << latticePosMax[2] << std::endl
              << "[smoothingFunctional] normalizer: "
              << normalizer << std::endl;
    return false;
  }

  // Normalizing to one
  for (int i=0; i<getSize(); i++) {
    _latticePosAndWeight[i].weight /= normalizer;
  }
  return true;
}


////////////////////// Class LinearAveragingSmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
LinearAveragingSmoothingFunctional<T, Lattice>::LinearAveragingSmoothingFunctional (
                 T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice )
         : SmoothingFunctional<T, Lattice>(kernelLength, converter, sLattice)
{}

template<typename T, typename Lattice>
T LinearAveragingSmoothingFunctional<T, Lattice>::compute(T physPosP[], T physPosL[])
{
  return this->smoothingFunction(physPosP[0] - physPosL[0])
       * this->smoothingFunction(physPosP[1] - physPosL[1])
       * this->smoothingFunction(physPosP[2] - physPosL[2]);
}


////////////////////// Class VolumeAveragingSmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
VolumeAveragingSmoothingFunctional<T, Lattice>::VolumeAveragingSmoothingFunctional (
                 T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice )
         : SmoothingFunctional<T, Lattice>(kernelLength, converter, sLattice)
{}

template<typename T, typename Lattice>
T VolumeAveragingSmoothingFunctional<T, Lattice>::compute(T physPosP[], T physPosL[])
{
  return this->smoothingFunction ( sqrt (
         pow(physPosP[0] - physPosL[0], 2) +
         pow(physPosP[1] - physPosL[1], 2) +
         pow(physPosP[2] - physPosL[2], 2) ) );
}


////////////////////// Class DeenSmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
DeenSmoothingFunctional<T, Lattice>::DeenSmoothingFunctional (
                 T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice )
         : LinearAveragingSmoothingFunctional<T, Lattice>(kernelLength, converter, sLattice)
{}

template<typename T, typename Lattice>
T DeenSmoothingFunctional<T, Lattice>::smoothingFunction(T delta)
{
  return ( pow(delta, 4)/pow(this->_kernelLength, 5)
           - 2.*pow(delta, 2)/pow(this->_kernelLength, 3)
           + 1./this->_kernelLength
         );
}


////////////////////// Class StepSmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
StepSmoothingFunctional<T, Lattice>::StepSmoothingFunctional (
                 T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice )
         : VolumeAveragingSmoothingFunctional<T, Lattice>(kernelLength, converter, sLattice)
{}

template<typename T, typename Lattice>
T StepSmoothingFunctional<T, Lattice>::smoothingFunction(T delta)
{
  return 1.;
}


}

#endif
