/*  DESCRIPTOR Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2016 Thomas Henn, Fabian Klemens, Robin Trunk, Davide Dapelo
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


#ifndef PARTICLEDYNAMICS_3D_HH
#define PARTICLEDYNAMICS_3D_HH

#include "hlbmDynamics3D.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::addCuboid(Vector< T, 3> center, T xLength, T yLength, T zLength, T density, T epsilon, Vector< T, 3 > theta, Vector<S,3> vel)
{
  _vectorOfIndicator.push_back (
		new SmoothIndicatorCuboid3D<T, T, true>(center, xLength, yLength, zLength, density, epsilon, theta, vel) );
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::addSphere(Vector< T, 3> center, T radius, T epsilon, T density, Vector<S,3> vel)
{
  _vectorOfIndicator.push_back (
		new SmoothIndicatorSphere3D<T, T, true>(center, radius, density, epsilon, vel) );
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::addParticle(SmoothIndicatorF3D<T, T, true>& indicator)
{
  _vectorOfIndicator.push_back(&indicator);
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::computeBoundaryForce(std::vector<SmoothIndicatorF3D<T,T,true>* >& indicator)
{
  SuperLatticePorousMomentumLossForce3D<T, DESCRIPTOR> force(_sLattice, _superGeometry, indicator, _converter);
  T sumF[force.getTargetDim()];
  for (int i=0; i<force.getTargetDim(); i++) {
    sumF[i]=0.;
  }

  int input[1];
  force(sumF, input);
  for (typename std::vector<SmoothIndicatorF3D<T, T, true>* >::size_type iInd=0; iInd!=indicator.size(); iInd++) {
    /// get particle acceleration through boundary force and gravity (and buoyancy)
    Vector<T,3> acceleration2;
    Vector<T,3> force;
    Vector<T,3> alpha2;
    force[0] = sumF[0+7*iInd];
    force[1] = sumF[1+7*iInd];
    force[2] = sumF[2+7*iInd];
    acceleration2[0] = sumF[0+7*iInd] / indicator[iInd]->getMass() + _accExt[0];
    acceleration2[1] = sumF[1+7*iInd] / indicator[iInd]->getMass() + _accExt[1];
    acceleration2[2] = sumF[2+7*iInd] / indicator[iInd]->getMass() + _accExt[2];
    alpha2[0] = sumF[3+7*iInd] / indicator[iInd]->getMofi()[0];
    alpha2[1] = sumF[4+7*iInd] / indicator[iInd]->getMofi()[1];
    alpha2[2] = sumF[5+7*iInd] / indicator[iInd]->getMofi()[2];
    indicator[iInd]->setForce( force );
    indicator[iInd]->setAcc2( acceleration2 );
    indicator[iInd]->setAlpha2( alpha2 );
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::addWallColl(SmoothIndicatorF3D<T, T, true>& indicator, T delta)
{
  T w1 = delta * .5;
  T w = delta * .5;

  T rad = indicator.getRadius();
  T massInv = 1. / indicator.getMass();

  std::vector<T> dx(3, T());
  dx[0] = _lengthX - _converter.getPhysDeltaX() - indicator.getPos()[0];
  dx[1] = _lengthY - _converter.getPhysDeltaX() - indicator.getPos()[1];
  dx[2] = indicator.getPos()[2] - _converter.getPhysDeltaX();

  for (int i = 0; i < 3; i++) {
    if (dx[i] <= rad) {
      indicator.getAcc2()[i] += massInv * -dx[i] * (rad - dx[i]) / w1;
    }
    if (indicator.getPos()[i] <= rad) {
      indicator.getAcc2()[i] += massInv * indicator.getPos()[i] * (rad - indicator.getPos()[i]) / w1;
    }
    if (dx[i] > rad && dx[i] <= rad + delta) {
      indicator.getAcc2()[i] += massInv * -dx[i] * std::pow((rad + delta - dx[i]), 2) / w;
    }
    if (indicator.getPos()[i] > rad && indicator.getPos()[i] <= rad + delta) {
      indicator.getAcc2()[i] += massInv * indicator.getPos()[i] * std::pow((rad + delta - indicator.getPos()[i]), 2) / w;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::verletIntegration(SmoothIndicatorF3D<T, T, true>& indicator)
{
  T time = _converter.getConversionFactorTime();
  T time2 = time*time;

  Vector<T,3> position, velocity, theta, omega, alpha;
  Vector<T,9> rotationMatrix; 
  for (int i=0; i<3; i++) {
    position[i] = indicator.getPos()[i] + indicator.getVel()[i] * time + (0.5 * indicator.getAcc()[i] * time2);
    T avgAcc = (indicator.getAcc()[i] +  indicator.getAcc2()[i]) * 0.5;
    velocity[i] = indicator.getVel()[i] + avgAcc * time;
    theta[i] = std::fmod( indicator.getTheta()[i] + indicator.getOmega()[i] * time + (0.5 * indicator.getAlpha()[i] * time2), 2.*M_PI );
    T avgAlpha = (indicator.getAlpha()[i] +  indicator.getAlpha2()[i]) * 0.5;
    omega[i] = indicator.getOmega()[i] + avgAlpha * time;
  }
  indicator.setPos( position );
  indicator.setVel( velocity );
  indicator.setAcc( indicator.getAcc2() );
  indicator.setTheta( theta );
  indicator.setOmega( omega );
  indicator.setAlpha( indicator.getAlpha2() );

  T cos0 = std::cos(indicator.getTheta()[0]);
  T cos1 = std::cos(indicator.getTheta()[1]);
  T cos2 = std::cos(indicator.getTheta()[2]);
  T sin0 = std::sin(indicator.getTheta()[0]);
  T sin1 = std::sin(indicator.getTheta()[1]);
  T sin2 = std::sin(indicator.getTheta()[2]);
  
  rotationMatrix[0] = cos1 * cos2;
  rotationMatrix[1] = sin0*sin1*cos2 - cos0*sin2;
  rotationMatrix[2] = cos0*sin1*cos2 + sin0*sin2;
  rotationMatrix[3] = cos1*sin2;
  rotationMatrix[4] = sin0*sin1*sin2 + cos0*cos2;
  rotationMatrix[5] = cos0*sin1*sin2 - sin0*cos2;
  rotationMatrix[6] = -sin1;
  rotationMatrix[7] = sin0*cos1;
  rotationMatrix[8] = cos0*cos1;
  indicator.setRotationMatrix( rotationMatrix );
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::updateParticleDynamics(std::string name, SmoothIndicatorF3D<T, T, true>& indicator)
{
  if (name == "euler") {
    this->eulerIntegration(indicator);
  } else if (name == "verlet") {
  this->verletIntegration(indicator);
  } else {
    std::cout << "ERROR: no valid integration...use 'euler' or 'verlet'"
              << std::endl;
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::checkAndRemoveEscaped()
{
  for (auto i=_vectorOfIndicator.begin(); i!=_vectorOfIndicator.end(); i++) {
    auto internal = std::make_shared<IndicatorLayer3D<T> > (*_indicatorF.get(), -_converter.getConversionFactorLength()+(**i).getEpsilon() );
    IndicMinus3D<T> layer(_indicatorF, internal);
    SuperLatticeIndicatorSmoothIndicatorIntersection3D<T,DESCRIPTOR,true> intersection(
      _sLattice, _superGeometry, layer, **i );
    int input[1];
    T output[1] = {0.};
    intersection(output, input);
    if (output[0] == 1) {
      _vectorOfIndicator.erase(i);
      if (i==_vectorOfIndicator.end() ) {
        break;
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::addParticleField(SmoothIndicatorF3D<T, T, true>& indicator)
{
  /// Analytical3D functor for particle motion (trans+rot)
  ParticleU3D<T,T,DESCRIPTOR> velocity(indicator, _converter);
  setSuperExternalParticleField(_superGeometry, velocity, indicator, _sLattice);
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::simulateTimestep(std::string name)
{
  // Remove particles from domain if _escapeFromDomain is toggled
  if (_escapeFromDomain) checkAndRemoveEscaped();

  // Compute force acting on particles boundary
  computeBoundaryForce(_vectorOfIndicator);

  // Update particle dynamics and porous particle field
  for (auto i=_vectorOfIndicator.begin(); i!=_vectorOfIndicator.end(); i++) {
    updateParticleDynamics(name, **i);
    addParticleField(**i);
  }
}


template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::print()
{
  OstreamManager clout(std::cout, "ParticleDynamics3D");
  clout << "Number of particles=" << _vectorOfIndicator.size() << std::endl;
  int count = 0;
  for (auto i=_vectorOfIndicator.begin(); i!=_vectorOfIndicator.end(); i++) {
    clout << "Particle " << ++count << " (" << (**i).name() << "):" << std::endl;
    clout << " |Circum radius(m)=     " << setw(13) << (**i).getCircumRadius() << std::endl;
    clout << " |Mass(kg)=             " << setw(13) << (**i).getMass() << std::endl;
    clout << " |Position(m)=         (" << setw(13) << (**i).getPos()[0] << ", " << setw(13) << (**i).getPos()[1] << ", " << setw(13) << (**i).getPos()[2] << ")" << std::endl;
    clout << " |Angle(°)=            (" << setw(13) << (**i).getTheta()[0]*(180/M_PI) << ", " << setw(13) << (**i).getTheta()[1]*(180/M_PI) << ", " << setw(13) << (**i).getTheta()[2]*(180/M_PI) << ")" << std::endl;
    clout << " |Velocity(m/s)=       (" << setw(13) << (**i).getVel()[0] << ", " << setw(13) << (**i).getVel()[1] << ", " << setw(13) << (**i).getVel()[2] << ")" << std::endl;
    clout << " |Ang. velocity(°/s)=  (" << setw(13) << (**i).getOmega()[0]*(180/M_PI) << ", " << setw(13) << (**i).getOmega()[1]*(180/M_PI) << ", " << setw(13) << (**i).getOmega()[2]*(180/M_PI) << ")" << std::endl;
    clout << " |Hydro. Force(N)=     (" << setw(13) << (**i).getForce()[0] << ", " << setw(13) << (**i).getForce()[1] << ", " << setw(13) << (**i).getForce()[2] << ")" << std::endl;
    clout << " |Acceleration(m/s^2)= (" << setw(13) << (**i).getAcc()[0] << ", " << setw(13) << (**i).getAcc()[1] << ", " << setw(13) << (**i).getAcc()[2] << ")" << std::endl;
    clout << " |Ang. acc.(°/s^2)=    (" << setw(13) << (**i).getAlpha()[0]*(180/M_PI) << ", " << setw(13) << (**i).getAlpha()[1]*(180/M_PI) << ", " << setw(13) << (**i).getAlpha()[1]*(180/M_PI) << ")" << std::endl;
  }
}

// WORKS ONLY FOR SPHERES FOR NOW!!
template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::load(std::string filename, T epsilon)
{
  std::ifstream iout( (singleton::directories().getLogOutDir() + filename).c_str() );

  while (iout) {
    std::string name = "";
    iout >> name;

    if (name == "sphere") {
      T radius = T(), mass = T();
      Vector<T, 3> center = {T(), T(), T()};
      Vector<T, 3> vel = {T(), T(), T()};
      iout >> center[0] >> center[1] >> center[2]
           >> radius >> mass
           >> vel[0] >> vel[1] >> vel[2];
      addSphere(center, radius, epsilon, mass, vel);
    }
  }

  iout.close();
}

// WORKS ONLY FOR SPHERES FOR NOW!!
template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::save(std::string filename)
{
  std::ofstream fout( (singleton::directories().getLogOutDir() + filename).c_str() );

  for (auto i=_vectorOfIndicator.begin(); i!=_vectorOfIndicator.end(); i++) {
    if ((**i).name() == "sphere") {
      fout << (**i).name()      << " "
           << (**i).getPos()[0] << " " << (**i).getPos()[1] << " " << (**i).getPos()[2] << " " 
           << (**i).getCircumRadius() << " " << (**i).getMass()   << " " 
           << (**i).getVel()[0] << " " << (**i).getVel()[1] << " " << (**i).getVel()[2]
           << std::endl;
    }
  }
  fout.close();
}

template<typename T, typename DESCRIPTOR>
void ParticleDynamics3D<T, DESCRIPTOR>::eulerIntegration(SmoothIndicatorF3D<T,T,true>& indicator) {
  T time = _converter.getConversionFactorTime();

  Vector<T,3> position, velocity, theta, omega, alpha;
  Vector<T,9> rotationMatrix; 
  for (int i=0; i<3; i++) {
    velocity[i] = indicator.getVel()[i] + indicator.getAcc2()[i] * time;
    position[i] = indicator.getPos()[i] + indicator.getVel()[i] * time;
    omega[i] = indicator.getOmega()[i] + indicator.getAlpha2()[i] * time;
    theta[i] = indicator.getTheta()[i] + indicator.getOmega()[i] * time;
  }

  indicator.setPos( position );
  indicator.setVel( velocity );
  indicator.setAcc( indicator.getAcc2() );
  indicator.setTheta( theta );
  indicator.setOmega( omega );
  indicator.setAlpha( indicator.getAlpha2() );

  T cos0 = std::cos(indicator.getTheta()[0]);
  T cos1 = std::cos(indicator.getTheta()[1]);
  T cos2 = std::cos(indicator.getTheta()[2]);
  T sin0 = std::sin(indicator.getTheta()[0]);
  T sin1 = std::sin(indicator.getTheta()[1]);
  T sin2 = std::sin(indicator.getTheta()[2]);
  
  rotationMatrix[0] = cos1 * cos2;
  rotationMatrix[1] = sin0*sin1*cos2 - cos0*sin2;
  rotationMatrix[2] = cos0*sin1*cos2 + sin0*sin2;
  rotationMatrix[3] = cos1*sin2;
  rotationMatrix[4] = sin0*sin1*sin2 + cos0*cos2;
  rotationMatrix[5] = cos0*sin1*sin2 - sin0*cos2;
  rotationMatrix[6] = -sin1;
  rotationMatrix[7] = sin0*cos1;
  rotationMatrix[8] = cos0*cos1;
  indicator.setRotationMatrix( rotationMatrix );
}

}

#endif
