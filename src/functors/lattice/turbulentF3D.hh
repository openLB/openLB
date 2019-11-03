/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Patrick Nathen, Mathias J. Krause
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

#ifndef TURBULENT_F_3D_HH
#define TURBULENT_F_3D_HH

#include<vector>
#include<cmath>

#include "turbulentF3D.h"
#include "blockLatticeLocalF3D.h"
#include "dynamics/smagorinskyBGKdynamics.h"
#include "core/superLattice3D.h"
#include "core/finiteDifference.h"
#include "geometry/superGeometry3D.h"
#include "utilities/vectorHelpers.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity


namespace olb {


///////////////////////////// SuperLatticeYplus3D //////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeYplus3D<T,DESCRIPTOR>::SuperLatticeYplus3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
    const UnitConverter<T,DESCRIPTOR>& converter, SuperGeometry3D<T>& superGeometry,
    IndicatorF3D<T>& indicator, const int material )
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,1),
    _superGeometry(superGeometry), _indicator(indicator), _material(material)
{
  this->getName() = "yPlus";
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeYplus3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];

  output[0]=T();
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<T> normalTemp(3,T());
    std::vector<T> normal(3,T());       // normalized
    T counter = T();
    T distance = T();
    if (_superGeometry.get(input) == 1) {
      for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
        if (_superGeometry.get(input[0],
                               input[1] + descriptors::c<DESCRIPTOR>(iPop,0),
                               input[2] + descriptors::c<DESCRIPTOR>(iPop,1),
                               input[3] + descriptors::c<DESCRIPTOR>(iPop,2)) == _material) {
          counter++;
          normalTemp[0] += descriptors::c<DESCRIPTOR>(iPop,0);
          normalTemp[1] += descriptors::c<DESCRIPTOR>(iPop,1);
          normalTemp[2] += descriptors::c<DESCRIPTOR>(iPop,2);
        }
      }
      if ( !util::nearZero(counter) ) {
        // get physical Coordinates at intersection

        std::vector<T> physR(3, T());
        _superGeometry.getCuboidGeometry().getPhysR(&(physR[0]), &(input[0]));

        T voxelSize = _superGeometry.getCuboidGeometry().get(globIC).getDeltaR();

        normal = util::normalize(normalTemp);

        std::vector<T> direction(normal);
        direction[0] = voxelSize*normal[0]*2.;
        direction[1] = voxelSize*normal[1]*2.;
        direction[2] = voxelSize*normal[2]*2.;

        // calculate distance to STL file
        if ( _indicator.distance(distance, physR, direction) ) {
          // call stress at this point
          T rho;
          T u[3];
          T pi[6];
          this->_sLattice.getBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix, lociy, lociz).computeRhoU(rho, u);
          this->_sLattice.getBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).computeStress(locix, lociy, lociz, pi);

          // Compute phys stress tau = mu*du/dx
          T omega = 1. / this->_converter.getLatticeRelaxationTime();
          T dt = this->_converter.getConversionFactorTime();
          T physFactor = -omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt*this->_converter.getPhysDensity(rho)*this->_converter.getPhysViscosity();

          //  Totel Stress projected from cell in normal direction on obstacle
          T Rx = pi[0]*physFactor*normal[0] + pi[1]*physFactor*normal[1] + pi[2]*physFactor*normal[2];
          T Ry = pi[1]*physFactor*normal[0] + pi[3]*physFactor*normal[1] + pi[4]*physFactor*normal[2];
          T Rz = pi[2]*physFactor*normal[0] + pi[4]*physFactor*normal[1] + pi[5]*physFactor*normal[2];

          // Stress appearing as pressure in corresponding direction is calculated and substracted
          T R_res_pressure = normal[0]*pi[0]*physFactor*normal[0] + normal[0]*pi[1]*physFactor*normal[1] + normal[0]*pi[2]*physFactor*normal[2]
                             +normal[1]*pi[1]*physFactor*normal[0] + normal[1]*pi[3]*physFactor*normal[1] + normal[1]*pi[4]*physFactor*normal[2]
                             +normal[2]*pi[2]*physFactor*normal[0] + normal[2]*pi[4]*physFactor*normal[1] + normal[2]*pi[5]*physFactor*normal[2];

          Rx -= R_res_pressure *normal[0];
          Ry -= R_res_pressure *normal[1];
          Rz -= R_res_pressure *normal[2];

          T tau_wall = sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
          T u_tau = sqrt(tau_wall/this->_converter.getPhysDensity(rho));
          //y_plus
          output[0] = distance*u_tau / this->_converter.getPhysViscosity();
        } // if 4
      }
    }
  }
  return true;
}

/*template <typename T, typename DESCRIPTOR>
BlockLatticeADM3D<T,DESCRIPTOR>::BlockLatticeADM3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, T sigma, int order, bool adaptive, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,5), _sigma(sigma), _order(order), _adaptive(adaptive), _converter(converter)
{
  this->getName() = "ADMfilter";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeADM3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  // Declaration of all variables needed for filtering

  int globX = input[0];
  int globY = input[1];
  int globZ = input[2];

  T output000[4] = {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ).computeRhoU(output000[0],output000+1);

  T outputp00[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX+1, globY, globZ).computeRhoU(outputp00[0],outputp00+1);

  T output2p00[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX+2, globY, globZ).computeRhoU(output2p00[0],output2p00+1);

  T outputn00[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX-1, globY, globZ).computeRhoU(outputn00[0],outputn00+1);

  T output2n00[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX-2, globY, globZ).computeRhoU(output2n00[0],output2n00+1);


  T output0p0[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY+1, globZ).computeRhoU(output0p0[0],output0p0+1);

  T output02p0[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY+2, globZ).computeRhoU(output02p0[0],output02p0+1);

  T output0n0[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY-1, globZ).computeRhoU(output0n0[0],output0n0+1);

  T output02n0[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY-2, globZ).computeRhoU(output02n0[0],output02n0+1);


  T output00p[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ+1).computeRhoU(output00p[0],output00p+1);

  T output002p[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ+2).computeRhoU(output002p[0],output002p+1);

  T output00n[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ-1).computeRhoU(output00n[0],output00n+1);

  T output002n[4]= {1.,0.,0.,0.};
  this-> _blockLattice.get(globX, globY, globZ-2).computeRhoU(output002n[0],output002n+1);

 T relax=_sigma;

   if (_adaptive==true ){
    ///////////////////////////////////////////////DISS VERSION///////////////////////////////////////////////////

  //  T diss = dissipation(u_000)[0];

    // std::cout << "diss: "<< diss << std::endl;

  //  T* avDiss = this-> _blockLattice.get(globX, globY , globZ)[localAvDissBeginsAt];

    // // std::cout <<"avDiss:" << *avDiss << std::endl;

  //   *avDiss = (*avDiss * this->_blockLattice.getStatistics().getTime() + diss) / (this->_blockLattice.getStatistics().getTime() + (int) 1 );

    // // std::cout <<"avDiss nach mittelung:" << *avDiss << std::endl;

   //  T TKE = 0.5*(velocity(u_000)[0]*velocity(u_000)[0]+velocity(u_000)[1]*velocity(u_000)[1]+velocity(u_000)[2]+velocity(u_000)[2]);

   //  T* avTKE = this-> _blockLattice.get(globX, globY , globZ)[localAvTKEBeginsAt];

   //  *avTKE = (*avTKE * this->_blockLattice.getStatistics().getTime() + TKE) / (this->_blockLattice.getStatistics().getTime() + (int) 1 );

    // std::cout << "TKE: "<< *avTKE << std::endl;


   //  relax = sqrt((diss - *avDiss)*(diss - *avDiss));// / (*avTKE * converter.getDeltaT());
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////// Stress Version ////////////////////////////////////////////////
    T stress[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
    this-> _blockLattice.computeStress(globX, globY, globZ, stress);

    T ux = stress(u_000)[0];
    T uy = stress(u_000)[1];
    T uz = stress(u_000)[2];

    T vx = stress(u_000)[3];
    T vy = stress(u_000)[4];
    T vz = stress(u_000)[5];

    T wx = stress(u_000)[6];
    T wy = stress(u_000)[7];
    T wz = stress(u_000)[8];


    T ux = stress[0];
    T uy = stress[1];
    T uz = stress[2];

    T vx = stress[3];
    T vy = stress[4];
    T vz = stress[5];

    T wx = stress[6];
    T wy = stress[7];
    T wz = stress[8];

    T norm = sqrt(  (wx*uz + vx*uy + ux*ux) + (wy*vz + vy*vy + uy*vx) + (wz*wz + vz*wy + uz*wx) ) ;


    T* avNorm = this-> _blockLattice.get(globX, globY , globZ)[_localAvDissBeginsAt];

    *avNorm = (*avNorm * this->_blockLattice.getStatistics().getTime() + norm) / (this->_blockLattice.getStatistics().getTime() + (int) 1 );

    // relax = sigma;
    // / (*avTKE * converter.getDeltaT());


    // if (this->_blockLattice.getStatistics().getTime() >= 30000){


     relax = sqrt((norm - *avNorm)*(norm - *avNorm)) ;

      // std::cout << "adaptive relaxation: " << relax <<  " time: "<<this->_blockLattice.getStatistics().getTime()<< endl;
    // }

   }

  if (_order==2) { // second order
    T d_0 = 6./16.;
    T d_1 = -4./16.;
    T d_2 = 1./16.;

    output[0] = output000[0] - relax*(d_2*(output2n00[0]+output02n0[0]+output002n[0]) +
                                       d_1*(outputn00[0]+output0n0[0]+output00n[0])+
                                       d_0*(output000[0]+output000[0]+output000[0])+
                                       d_1*(outputp00[0]+output0p0[0]+output00p[0])+
                                       d_2*(output2p00[0]+output02p0[0]+output002p[0]) );
    for (int i = 1; i < 4; ++i ) {
      output[i] = output000[i] - relax*(d_2*(output2n00[i] + output02n0[i] + output002n[i]) +
                                         d_1*(outputn00[i] + output0n0[i] + output00n[i])+
                                         d_0*(output000[i] + output000[i] + output000[i])+
                                         d_1*(outputp00[i] + output0p0[i] + output00p[i])+
                                         d_2*(output2p00[i] + output02p0[i] + output002p[i]) );
    }
  } else { // third order

    T output3p00[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX+3, globY, globZ).computeRhoU(output3p00[0],output3p00+1);

    T output3n00[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX-3, globY, globZ).computeRhoU(output3n00[0],output3n00+1);

    T output03p0[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX, globY+3, globZ).computeRhoU(output03p0[0],output03p0+1);

    T output03n0[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX, globY-3, globZ).computeRhoU(output03n0[0],output03n0+1);

    T output003p[4]= {1.,0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ+3).computeRhoU(output003p[0],output003p+1);

    T output003n[4] = {1.,0.,0.,0.};
    this-> _blockLattice.get(globX, globY, globZ-3).computeRhoU(output003n[0],output003n+1);

    T d_0 = 5./16.;
    T d_1 = -15./64.;
    T d_2 = 3./32.;
    T d_3 = -1./64.;
    output[0] = output000[0] - _sigma*(d_3*(output3n00[0]+output03n0[0]+output003n[0])+
                                       d_2*(output2n00[0]+output02n0[0]+output002n[0]) +
                                       d_1*(outputn00[0]+output0n0[0]+output00n[0])+
                                       d_0*(output000[0]+output000[0]+output000[0])+
                                       d_1*(outputp00[0]+output0p0[0]+output00p[0])+
                                       d_2*(output2p00[0]+output02p0[0]+output002p[0])+
                                       d_3*(output3p00[0]+output03p0[0]+output003p[0]) );
    for (int i = 1; i < 4; ++i ) {
      output[i] = output000[i] - _sigma*(d_3*(output3n00[i]+output03n0[i]+output003n[i])+
                                         d_2*(output2n00[i]+output02n0[i]+output002n[i]) +
                                         d_1*(outputn00[i]+output0n0[i]+output00n[i])+
                                         d_0*(output000[i]+output000[i]+output000[i])+
                                         d_1*(outputp00[i]+output0p0[i]+output00p[i])+
                                         d_2*(output2p00[i]+output02p0[i]+output002p[i])+
                                         d_3*(output3p00[i]+output03p0[i]+output003p[i]) );
    }
  }
  output[4]=relax;
  return true;
}

template <typename T, typename DESCRIPTOR>
void BlockLatticeADM3D<T,DESCRIPTOR>::execute(const int input[])
{
  T output[5] = {T(),T(),T(),T(),T()};
  this->operator()(output,input);
  this-> _blockLattice.get(input[0],input[1],input[2]).defineField<descriptors::FIL_RHO>( &output[0] );
  this-> _blockLattice.get(input[0],input[1],input[2]).defineField<descriptors::LOCAL_FIL_VEL_X>( &output[1]);
  this-> _blockLattice.get(input[0],input[1],input[2]).defineField<descriptors::LOCAL_FIL_VEL_Y>( &output[2]);
  this-> _blockLattice.get(input[0],input[1],input[2]).defineField<descriptors::LOCAL_FIL_VEL_Z>( &output[3]);
  this-> _blockLattice.get(input[0],input[1],input[2]).defineField<descriptors::LOCAL_SIGMA_ADM>( &output[4]);
}

template <typename T, typename DESCRIPTOR>
void BlockLatticeADM3D<T,DESCRIPTOR>::execute()
{
  int nX = this-> _blockLattice.getNx();
  int nY = this-> _blockLattice.getNy();
  int nZ = this-> _blockLattice.getNz();
  int i[3];
  for (i[0]=0; i[0]<nX; ++i[0]) {
    for (i[1]=0; i[1]<nY; ++i[1]) {
      for (i[2]=0; i[2]<nZ; ++i[2]) {
        execute(i);
      }
    }
  }
}


template <typename T, typename DESCRIPTOR>
SuperLatticeADM3D<T,DESCRIPTOR>::SuperLatticeADM3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, T sigma, int order, bool adaptive, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,5), _sigma(sigma), _order(order), _adaptive(adaptive), _converter(converter)
{
  this->getName() = "ADMfilter";
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeADM3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    int inputLocal[3]= {};
    T overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticeADM3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), _sigma, _order,_adaptive, _converter );

    blockLatticeF(output,inputLocal);
    return true;
  } else {
    return false;
  }
}

template <typename T, typename DESCRIPTOR>
void SuperLatticeADM3D<T,DESCRIPTOR>::execute(SuperGeometry3D<T>& superGeometry, const int material)
{
  this->_sLattice.communicate();
  CuboidGeometry3D<T>& cGeometry =  this->_sLattice.getCuboidGeometry();
  LoadBalancer<T>& load = this->_sLattice.getLoadBalancer();
  int overlap = this->_sLattice.getOverlap();

  for (int iC = 0; iC < load.size(); ++iC) {
    BlockLatticeADM3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(iC), _sigma, _order, _adaptive, _converter );
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    int nZ = cGeometry.get(load.glob(iC)).getNz();

    int i[3];
    for (i[0]=overlap; i[0]<nX+overlap; ++i[0]) {
      for (i[1]=overlap; i[1]<nY+overlap; ++i[1]) {
        for (i[2]=overlap; i[2]<nZ+overlap; ++i[2]) {
          if (superGeometry.getExtendedBlockGeometry(iC).get(i[0],i[1],i[2]) == material) {
            blockLatticeF.execute(i);
          }
        }
      }
    }
  }
}

*/

////////////////////////BlockFiniteDifference3D//////////////////////////////////
template <typename T>
BlockFiniteDifference3D<T>::BlockFiniteDifference3D
(BlockGeometryStructure3D<T>& blockGeometry, BlockF3D<T>& blockFunctor, std::list<int>& matNumber)
  : BlockF3D<T>(blockFunctor.getBlockStructure(), 3*blockFunctor.getTargetDim()), _blockGeometry(blockGeometry), _blockFunctor(blockFunctor), _matNumber(matNumber)
{
  this->getName() = "FiniteDifference";
  _targetDim = _blockFunctor.getTargetDim();
  _n[0] = this-> _blockGeometry.getNx()-1;
  _n[1] = this-> _blockGeometry.getNy()-1;
  _n[2] = this-> _blockGeometry.getNz()-1;

}

template <typename T>
bool BlockFiniteDifference3D<T>::operator() (T output[], const int input[])
{
//  // derivation tensor
  std::vector<std::vector<T>> fdGrad;

  fdGrad.resize(_targetDim);
  for (int i = 0; i < _targetDim; i++) {
    fdGrad[i].resize(3);
  }

  for (int i = 0; i < 3; i++) {
    int fInput_p[3];
    fInput_p[0] = input[0];
    fInput_p[1] = input[1];
    fInput_p[2] = input[2];
    fInput_p[i]+=1;

    int fInput_2p[3];
    fInput_2p[0] = input[0];
    fInput_2p[1] = input[1];
    fInput_2p[2] = input[2];
    fInput_2p[i]+=2;

    int fInput_3p[3];
    fInput_3p[0] = input[0];
    fInput_3p[1] = input[1];
    fInput_3p[2] = input[2];
    fInput_3p[i]+=3;

    int fInput_4p[3];
    fInput_4p[0] = input[0];
    fInput_4p[1] = input[1];
    fInput_4p[2] = input[2];
    fInput_4p[i]+=4;

    int fInput_n[3];
    fInput_n[0] = input[0];
    fInput_n[1] = input[1];
    fInput_n[2] = input[2];
    fInput_n[i]-=1;

    int fInput_2n[3];
    fInput_2n[0] = input[0];
    fInput_2n[1] = input[1];
    fInput_2n[2] = input[2];
    fInput_2n[i]-=2;

    int fInput_3n[3];
    fInput_3n[0] = input[0];
    fInput_3n[1] = input[1];
    fInput_3n[2] = input[2];
    fInput_3n[i]-=3;

    int fInput_4n[3];
    fInput_4n[0] = input[0];
    fInput_4n[1] = input[1];
    fInput_4n[2] = input[2];
    fInput_4n[i]-=4;

    T fOutput[_targetDim];
    _blockFunctor(fOutput,input);

    if (input[i] < 3) {
      if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get(fInput_2p[0], fInput_2p[1], fInput_2p[2])) == _matNumber.end()) {
        T fOutput_p[_targetDim];
        _blockFunctor(fOutput_p,fInput_p);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]= -fOutput[j] + fOutput_p[j];
        }
      } else {
        T fOutput_p[_targetDim];
        _blockFunctor(fOutput_p,fInput_p);
        T fOutput_2p[_targetDim];
        _blockFunctor(fOutput_2p,fInput_2p);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]=fd::boundaryGradient(fOutput[j], fOutput_p[j], fOutput_2p[j]);
        }
      }
    } else if (input[i] > _n[i]-3) {
      if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get(fInput_2n[0], fInput_2n[1], fInput_2n[2])) == _matNumber.end()) {
        T fOutput_n[_targetDim];
        _blockFunctor(fOutput_n,fInput_n);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]= -fOutput_n[j] + fOutput[j];
        }
      } else {
        T fOutput_n[_targetDim];
        _blockFunctor(fOutput_n,fInput_n);
        T fOutput_2n[_targetDim];
        _blockFunctor(fOutput_2n,fInput_2n);
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]=fd::boundaryGradient(-fOutput[j], -fOutput_n[j], -fOutput_2n[j]);
        }
      }
    } else {
      if ( std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get(fInput_n[0], fInput_n[1], fInput_n[2])) == _matNumber.end()  &&
           std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get(fInput_p[0], fInput_p[1], fInput_p[2])) == _matNumber.end() ) {
        for (int j=0; j < _targetDim; j++) {
          fdGrad[j][i]=0.;
        }
        // boundary treatment with Second-order asymmetric gradient
      } else if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get(fInput_n[0], fInput_n[1], fInput_n[2])) == _matNumber.end()) {
        if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get(fInput_2p[0], fInput_2p[1], fInput_2p[2])) == _matNumber.end()) {
          T fOutput_p[_targetDim];
          _blockFunctor(fOutput_p,fInput_p);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]= -fOutput[j] + fOutput_p[j];
          }
        } else {
          T fOutput_p[_targetDim];
          _blockFunctor(fOutput_p,fInput_p);
          T fOutput_2p[_targetDim];
          _blockFunctor(fOutput_2p,fInput_2p);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]=fd::boundaryGradient(fOutput[j], fOutput_p[j], fOutput_2p[j]);
          }
        }
      } else if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get(fInput_p[0], fInput_p[1], fInput_p[2])) == _matNumber.end() ) {
        if (std::find(_matNumber.begin(), _matNumber.end(), _blockGeometry.get(fInput_2n[0], fInput_2n[1], fInput_2n[2])) == _matNumber.end()) {
          T fOutput_n[_targetDim];
          _blockFunctor(fOutput_n,fInput_n);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]= -fOutput_n[j] + fOutput[j];
          }
        } else {
          T fOutput_n[_targetDim];
          _blockFunctor(fOutput_n,fInput_n);
          T fOutput_2n[_targetDim];
          _blockFunctor(fOutput_2n,fInput_2n);
          for (int j=0; j < _targetDim; j++) {
            fdGrad[j][i]=fd::boundaryGradient(-fOutput[j], -fOutput_n[j], -fOutput_2n[j]);
          }
        }
      } else {
        //inner domain 8th order central difference
        T fOutput_n[_targetDim];
        _blockFunctor(fOutput_n,fInput_n);

        T fOutput_2n[_targetDim];
        _blockFunctor(fOutput_2n,fInput_2n);

        T fOutput_3n[_targetDim];
        _blockFunctor(fOutput_3n,fInput_3n);

        T fOutput_4n[_targetDim];
        _blockFunctor(fOutput_4n,fInput_4n);

        T fOutput_p[_targetDim];
        _blockFunctor(fOutput_p,fInput_p);

        T fOutput_2p[_targetDim];
        _blockFunctor(fOutput_2p,fInput_2p);

        T fOutput_3p[_targetDim];
        _blockFunctor(fOutput_3p,fInput_3p);

        T fOutput_4p[_targetDim];
        _blockFunctor(fOutput_4p,fInput_4p);
        for (int j=0; j < _targetDim; j++) {
          //fdGrad[j][i]=fd::centralGradient(fOutput_p[j], fOutput_n[j]);
          fdGrad[j][i]=((T)672*(fOutput_p[j]-fOutput_n[j])+(T)168*(fOutput_2n[j]-fOutput_2p[j])
                        +(T)32*(fOutput_3p[j]-fOutput_3n[j])+(T)3*(fOutput_4n[j]-fOutput_4p[j])) / 840.;
        }
      }
    }
    for (int i=0; i < 3; i++) {
      for (int j=0; j < _targetDim; j++) {
        output[i*3+j] = fdGrad[i][j];
      }
    }
  }
  return true;
}



////////////////////////SuperFiniteDifference3D//////////////////////////////////
template <typename T>
SuperFiniteDifference3D<T>::SuperFiniteDifference3D
(SuperGeometry3D<T>& sGeometry, SuperF3D<T>& sFunctor, std::list<int>& matNumber) : SuperF3D<T>(sFunctor.getSuperStructure(),3*sFunctor.getTargetDim()),
  _sGeometry(sGeometry) ,_sFunctor(sFunctor), _matNumber(matNumber)
{
  this->getName() = "FiniteDifference";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockFiniteDifference3D<T> ( _sGeometry.getBlockGeometry(iC), _sFunctor.getBlockF(iC), _matNumber ));
  }
}

template <typename T>
bool SuperFiniteDifference3D<T>::operator() (T output[], const int input[])
{
  if (this->_superStructure.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_superStructure.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

////////////////////////BlockPhysFiniteDifference3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockPhysFiniteDifference3D<T,DESCRIPTOR>::BlockPhysFiniteDifference3D
(BlockF3D<T>& blockFinDiff, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockF3D<T>(blockFinDiff.getBlockStructure(), 3*blockFinDiff.getTargetDim()), _blockFinDiff(blockFinDiff), _converter(converter)
{
  this->getName() = "PhysFiniteDifference";
  _targetDim = _blockFinDiff.getTargetDim();

}

template <typename T, typename DESCRIPTOR>
bool BlockPhysFiniteDifference3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  _blockFinDiff(output,input);
  for (int i = 0; i < _targetDim; i++) {
    output[i] /= _converter.getConversionFactorLength();
  }
  return true;
}



////////////////////////SuperPhysFiniteDifference3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperPhysFiniteDifference3D<T,DESCRIPTOR>::SuperPhysFiniteDifference3D
(SuperGeometry3D<T>& sGeometry, SuperF3D<T>& sFunctor, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperF3D<T>(sFunctor.getSuperStructure(),3*sFunctor.getTargetDim()),
  _sFinDiff(sGeometry,sFunctor,matNumber),_converter(converter)
{
  this->getName() = "PhysFiniteDifference";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockPhysFiniteDifference3D<T,DESCRIPTOR> (_sFinDiff.getBlockF(iC), _converter ));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperPhysFiniteDifference3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_superStructure.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_superStructure.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}
////////////////////////BlockLatticeVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeVelocityGradientFD3D<T,DESCRIPTOR>::BlockLatticeVelocityGradientFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFinDiff)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockFinDiff(blockFinDiff)
{
  this->getName() = "VelocityGradientFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeVelocityGradientFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //1 dudx 2 dudy 3 dudz
  //4 dydx 5 dydy 6 dydz
  //7 dwdx 8 dwdy 9 dwdz
  _blockFinDiff(output,input);
  return true;
}

////////////////////////BlockLatticeExternalVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeExternalVelocityGradientFD3D<T,DESCRIPTOR>::BlockLatticeExternalVelocityGradientFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFinDiff)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockFinDiff(blockFinDiff)
{
  this->getName() = "externalVelocityGradientFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeExternalVelocityGradientFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //1 dudx 2 dudy 3 dudz
  //4 dydx 5 dydy 6 dydz
  //7 dwdx 8 dwdy 9 dwdz
  _blockFinDiff(output,input);
  return true;
}

////////////////////////SuperLatticeVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeVelocityGradientFD3D<T,DESCRIPTOR>::SuperLatticeVelocityGradientFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVelocity(sLattice), _sFinDiff(sGeometry, _sVelocity, matNumber)
{
  this->getName() = "VelocityGradientFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVelocityGradientFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sFinDiff.getBlockF(iC)));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeVelocityGradientFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}
////////////////////////SuperLatticeExternalVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeExternalVelocityGradientFD3D<T,DESCRIPTOR>::SuperLatticeExternalVelocityGradientFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVelocity(sLattice), _sFinDiff(sGeometry, _sVelocity, matNumber)
{
  this->getName() = "externalVelocityGradientFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeExternalVelocityGradientFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sFinDiff.getBlockF(iC)));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeExternalVelocityGradientFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}
////////////////////////BlockLatticePhysVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysVelocityGradientFD3D<T,DESCRIPTOR>::BlockLatticePhysVelocityGradientFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFinDiff, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockFinDiff(blockFinDiff), _converter(converter)
{
  this->getName() = "PhysVelocityGradientFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysVelocityGradientFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  _blockFinDiff(output,input);
  return true;
}



////////////////////////SuperLatticePhysVelocityGradientFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysVelocityGradientFD3D<T,DESCRIPTOR>::SuperLatticePhysVelocityGradientFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVelocity(sLattice, converter), _sFinDiff(sGeometry, _sVelocity, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysVelocityGradientFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysVelocityGradientFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sFinDiff.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePhysVelocityGradientFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}
////////////////////////BlockLatticeStrainRateFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeStrainRateFD3D<T,DESCRIPTOR>::BlockLatticeStrainRateFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockVeloGrad(blockVeloGrad)
{
  this->getName() = "StrainRateFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeStrainRateFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0];
  output[1] = 0.5 * (velograd[1] + velograd[3]);
  output[2] = 0.5 * (velograd[2] + velograd[6]);
  output[3] = output[1];
  output[4] = velograd[4];
  output[5] = 0.5 * (velograd[5] + velograd[7]);
  output[6] = output[2];
  output[7] = output[5];
  output[8] = velograd[8];
  return true;
}



////////////////////////SuperLatticeStrainRateFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeStrainRateFD3D<T,DESCRIPTOR>::SuperLatticeStrainRateFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVeloGrad(sGeometry, sLattice, matNumber)
{
  this->getName() = "StrainRateFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeStrainRateFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeStrainRateFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}
////////////////////////BlockLatticePhysStrainRateFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysStrainRateFD3D<T,DESCRIPTOR>::BlockLatticePhysStrainRateFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "PhysStrainRateFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysStrainRateFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0];
  output[1] = 0.5 * (velograd[1] + velograd[3]);
  output[2] = 0.5 * (velograd[2] + velograd[6]);
  output[3] = output[1];
  output[4] = velograd[4];
  output[5] = 0.5 * (velograd[5] + velograd[7]);
  output[6] = output[2];
  output[7] = output[5];
  output[8] = velograd[8];

  return true;
}



////////////////////////SuperLatticePhysStrainRateFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysStrainRateFD3D<T,DESCRIPTOR>::SuperLatticePhysStrainRateFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysStrainRateFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysStrainRateFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePhysStrainRateFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

////////////////////////BlockLatticeDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeDissipationFD3D<T,DESCRIPTOR>::BlockLatticeDissipationFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "DissipationFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0] * velograd[0] + velograd[1] * velograd[1] + velograd[2] * velograd[2] +
              velograd[3] * velograd[3] + velograd[4] * velograd[4] + velograd[5] * velograd[5] +
              velograd[6] * velograd[6] + velograd[7] * velograd[7] + velograd[8] * velograd[8];
  output[0] *= _converter.getLatticeViscosity();

  return true;
}



////////////////////////SuperLatticeDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeDissipationFD3D<T,DESCRIPTOR>::SuperLatticeDissipationFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "DissipationFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeDissipationFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

////////////////////////BlockLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysDissipationFD3D<T,DESCRIPTOR>::BlockLatticePhysDissipationFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "PhysDissipationFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0] * velograd[0] + velograd[1] * velograd[1] + velograd[2] * velograd[2] +
              velograd[3] * velograd[3] + velograd[4] * velograd[4] + velograd[5] * velograd[5] +
              velograd[6] * velograd[6] + velograd[7] * velograd[7] + velograd[8] * velograd[8];
  output[0] *= _converter.getPhysViscosity();

  return true;
}



////////////////////////SuperLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysDissipationFD3D<T,DESCRIPTOR>::SuperLatticePhysDissipationFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysDissipationFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysDissipationFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePhysDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

////////////////////////BlockLatticeEffectiveDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeEffectiveDissipationFD3D<T,DESCRIPTOR>::BlockLatticeEffectiveDissipationFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVeloGrad(blockVeloGrad), _converter(converter), _LESdynamics(LESdynamics)
{
  this->getName() = "EffectiveDissipationFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeEffectiveDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0] * velograd[0] + velograd[1] * velograd[1] + velograd[2] * velograd[2] +
              velograd[3] * velograd[3] + velograd[4] * velograd[4] + velograd[5] * velograd[5] +
              velograd[6] * velograd[6] + velograd[7] * velograd[7] + velograd[8] * velograd[8];

  T omegaEff = _LESdynamics.getEffectiveOmega(this->_blockLattice.get(input[0], input[1], input[2]));
  T nuEff = ((1./omegaEff)-0.5)/descriptors::invCs2<T,DESCRIPTOR>();
  output[0] *= nuEff;

  return true;
}



////////////////////////SuperLatticeEffectiveDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeEffectiveDissipationFD3D<T,DESCRIPTOR>::SuperLatticeEffectiveDissipationFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1), _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "EffectiveDissipationFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeEffectiveDissipationFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeEffectiveDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

////////////////////////BlockLatticePhysEffectiveDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::BlockLatticePhysEffectiveDissipationFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVeloGrad(blockVeloGrad), _converter(converter), _LESdynamics(LESdynamics)
{
  this->getName() = "PhysEffectiveDissipationFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[0] * velograd[0] + velograd[1] * velograd[1] + velograd[2] * velograd[2] +
              velograd[3] * velograd[3] + velograd[4] * velograd[4] + velograd[5] * velograd[5] +
              velograd[6] * velograd[6] + velograd[7] * velograd[7] + velograd[8] * velograd[8];

  T omegaEff = _LESdynamics.getEffectiveOmega(this->_blockLattice.get(input[0], input[1], input[2]));
  T nuEff = ((1./omegaEff)-0.5)/descriptors::invCs2<T,DESCRIPTOR>();
  output[0] *= _converter.getPhysViscosity( nuEff );

  return true;
}



////////////////////////SuperLatticePhysEffectiveDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::SuperLatticePhysEffectiveDissipationFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1), _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysEffectiveDissipationFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sVeloGrad.getBlockF(iC), this->_converter, LESdynamics));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePhysEffectiveDissipationFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

////////////////////////BlockLatticeVorticityFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeVorticityFD3D<T,DESCRIPTOR>::BlockLatticeVorticityFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 3), _blockVeloGrad(blockVeloGrad)
{
  this->getName() = "VorticityFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeVorticityFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[7] - velograd[5];
  output[1] = velograd[2] - velograd[6];
  output[2] = velograd[3] - velograd[1];
  return true;
}



////////////////////////SuperLatticeVorticityFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeVorticityFD3D<T,DESCRIPTOR>::SuperLatticeVorticityFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3),
  _sVeloGrad(sGeometry, sLattice, matNumber)
{
  this->getName() = "VorticityFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVorticityFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeVorticityFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

////////////////////////BlockLatticePhysVorticityFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysVorticityFD3D<T,DESCRIPTOR>::BlockLatticePhysVorticityFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 3), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "PhysVorticityFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysVorticityFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
  output[0] = velograd[7] - velograd[5];
  output[1] = velograd[2] - velograd[6];
  output[2] = velograd[3] - velograd[1];
  return true;
}



////////////////////////SuperLatticePhysVorticityFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysVorticityFD3D<T,DESCRIPTOR>::SuperLatticePhysVorticityFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysVorticityFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysVorticityFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePhysVorticityFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

////////////////////////BlockLatticePhysStressFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysStressFD3D<T,DESCRIPTOR>::BlockLatticePhysStressFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockStrainRate, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 9), _blockStrainRate(blockStrainRate), _converter(converter)
{
  this->getName() = "PhysStressFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysStressFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  _blockStrainRate(output,input);
  for (int i = 0; i < 9; i++) {
    output[i] /= _converter.getPhysViscosity() * _converter.getPhysDensity();
  }
  return true;
}



////////////////////////SuperLatticePhysStressFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysStressFD3D<T,DESCRIPTOR>::SuperLatticePhysStressFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9),
  _sStrainRate(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysStressFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysStressFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sStrainRate.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePhysStressFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

////////////////////////BlockIsotropicHomogeneousTKE//////////////////////////////////
template<typename T, typename DESCRIPTOR>
BlockIsotropicHomogeneousTKE3D<T, DESCRIPTOR>::BlockIsotropicHomogeneousTKE3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVelocity)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _blockVelocity(blockVelocity)
{
  this->getName() = "IsotropicHomogeneousTKE";
}

template<typename T, typename DESCRIPTOR>
bool BlockIsotropicHomogeneousTKE3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();
  T data[_blockVelocity.getTargetDim()];
  _blockVelocity(data,input);
  for (int i = 0; i < _blockVelocity.getTargetDim(); ++i) {
    output[0] +=  data[i] * data[i];
  }
  output[0] = 0.5 * output[0];
  return true;
}


////////////////////////SuperIsotropicHomogeneousTKE//////////////////////////////////
template<typename T, typename DESCRIPTOR>
SuperIsotropicHomogeneousTKE3D<T, DESCRIPTOR>::SuperIsotropicHomogeneousTKE3D(
SuperLattice3D<T,DESCRIPTOR>& sLattice, const  UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, 1), _sVelocity(sLattice, converter), _converter(converter)

{
  this->getName() = "IsotropicHomogeneousTKE";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++)
  {
    this->_blockF.emplace_back(new BlockIsotropicHomogeneousTKE3D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this-> _sVelocity.getBlockF(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperIsotropicHomogeneousTKE3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

 ////////////////////////BlockLatticePhysEnstrophyFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticePhysEnstrophyFD3D<T,DESCRIPTOR>::BlockLatticePhysEnstrophyFD3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const  UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 3), _blockVeloGrad(blockVeloGrad), _converter(converter)
{
  this->getName() = "PhysEnstrophyFD";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysEnstrophyFD3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T velograd[9];
  _blockVeloGrad(velograd,input);
    // mod 2019-04-10:
    // this is now the local enstrophy value: 0.5 * vort_\alpha * vort_\alpha.
    // to obtain global enstrophy: integrate and divide by integration volume.
    output[0] = 0.5 * ( pow(velograd[7] - velograd[5], 2) + pow(velograd[2] - velograd[6], 2) + pow(velograd[3] - velograd[1], 2) );
    // output[1] = pow(velograd[2] - velograd[6], 2);
    // output[2] = pow(velograd[3] - velograd[1], 2);
  return true;
}

////////////////////////SuperLatticePhysEnstrophyFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticePhysEnstrophyFD3D<T,DESCRIPTOR>::SuperLatticePhysEnstrophyFD3D
(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const  UnitConverter<T,DESCRIPTOR>& converter) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3),
  _sVeloGrad(sGeometry, sLattice, matNumber, converter), _converter(converter)
{
  this->getName() = "PhysEnstrophyFD";
  int maxC = this->_superStructure.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysEnstrophyFD3D<T,DESCRIPTOR> (this->_sLattice.getBlockLattice(iC), this->_sVeloGrad.getBlockF(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePhysEnstrophyFD3D<T, DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

/*
////////////////////////BlockLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
BlockLatticeSigmaADM3D<T,DESCRIPTOR>::BlockLatticeSigmaADM3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice )
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,3)
{
  this->getName() = "SigmaADM3D";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeSigmaADM3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

 T* sigma = this-> _blockLattice.get(input[0], input[1] , input[2])[_localSigmaADMBeginsAt];
   output[0] = *sigma;

  return true;

}


////////////////////////SuperLatticePhysDissipationFD3D//////////////////////////////////
template <typename T, typename DESCRIPTOR>
SuperLatticeSigmaADM3D<T,DESCRIPTOR>::SuperLatticeSigmaADM3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3)
{
  this->getName() = "SigmaADM3D";
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeSigmaADM3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    int inputLocal[3]= {};
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticeSigmaADM3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)));

    blockLatticeF(output,inputLocal);
    return true;
  } else {
    return false;
  }


}

*/

} // end namespace olb
#endif
