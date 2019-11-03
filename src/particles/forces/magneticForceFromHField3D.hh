/*
 *  Copyright (C) 2018 Marie-Luise Maier
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

#ifndef MagneticForceFromHField3D_HH
#define MagneticForceFromHField3D_HH

#include <cmath>
#include <vector>
#include "magneticForceFromHField3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

/*
template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
MagneticForceFromHField3D<T, PARTICLETYPE, DESCRIPTOR>::MagneticForceFromHField3D(
    SuperLattice3D<T, DESCRIPTOR>& sLattice,
    SuperLatticeField3D<T,DESCRIPTOR>& sLatticeHField, T Mp, T scale) :
    Force3D<T, PARTICLETYPE>(), _sLattice(sLattice), _sLatticeHField(
        sLatticeHField), _Mp(Mp), _scale(scale) {
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void MagneticForceFromHField3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
    typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
    ParticleSystem3D<T, PARTICLETYPE>& pSys) {

  // vacuum permeability
  T mu_0 = 4. * M_PI * 1.e-7;
  T Vp = 4. / 3. * M_PI * std::pow(p->getRad(), 3);
  T physLatticeLength =
      _sLattice.getCuboidGeometry().get(p->getCuboid()).getDeltaR();

  int glob = p->getCuboid();
  int latticePos[3] = { 0, 0, 0 };
  _sLattice.getCuboidGeometry().get(p->getCuboid()).getLatticeR(latticePos,
      &p->getPos()[0]);

  T H[3] = { T(), T(), T() };
  T HXplus[3] = { T(), T(), T() };
  T HYplus[3] = { T(), T(), T() };
  T HZplus[3] = { T(), T(), T() };

  int X[4] = { glob, latticePos[0], latticePos[1], latticePos[2] };

  int Xplus[4] = { glob, latticePos[0] + 1, latticePos[1], latticePos[2] };
  int Yplus[4] = { glob, latticePos[0], latticePos[1] + 1, latticePos[2] };
  int Zplus[4] = { glob, latticePos[0], latticePos[1], latticePos[2] + 1 };

  T gradHx[3] = { T(), T(), T() };
  T gradHy[3] = { T(), T(), T() };
  T gradHz[3] = { T(), T(), T() };
  T modGradH[3] = { T(), T(), T() };

  _sLatticeHField(H, X);
  _sLatticeHField(HXplus, Xplus);
  _sLatticeHField(HYplus, Yplus);
  _sLatticeHField(HZplus, Zplus);

  if (!(util::nearZero(H[0]) && util::nearZero(H[1]) && util::nearZero(H[2]))) {
    Vector<T, 3> M(H[0], H[1], H[2]);
    M.normalize();
    M[0] *= _Mp;
    M[1] *= _Mp;
    M[2] *= _Mp;

    gradHx[0] = (HXplus[0] - H[0]) / physLatticeLength;
    gradHx[1] = (HXplus[1] - H[1]) / physLatticeLength;
    gradHx[2] = (HXplus[2] - H[2]) / physLatticeLength;

    gradHy[0] = (HYplus[0] - H[0]) / physLatticeLength;
    gradHy[1] = (HYplus[1] - H[1]) / physLatticeLength;
    gradHy[2] = (HYplus[2] - H[2]) / physLatticeLength;

    gradHz[0] = (HZplus[0] - H[0]) / physLatticeLength;
    gradHz[1] = (HZplus[1] - H[1]) / physLatticeLength;
    gradHz[2] = (HZplus[2] - H[2]) / physLatticeLength;

    // modGradH = (M \cdot \nabla) H
    modGradH[0] = M[0] * gradHx[0] + M[1] * gradHy[0] + M[2] * gradHz[0];
    modGradH[1] = M[0] * gradHx[1] + M[1] * gradHy[1] + M[2] * gradHz[1];
    modGradH[2] = M[0] * gradHx[2] + M[1] * gradHy[2] + M[2] * gradHz[2];

    if (modGradH[0] > T() || modGradH[1] > T() || modGradH[2] > T()) {
      std::cout << "without scale: modGradH[0] " << Vp * mu_0 * modGradH[0]
          << "modGradH[1] " << Vp * mu_0 * modGradH[1] << "modGradH[2] "
          << Vp * mu_0 * modGradH[2] << std::endl;
    }
    p->getForce()[0] += Vp * mu_0 * modGradH[0] * _scale;
    p->getForce()[1] += Vp * mu_0 * modGradH[1] * _scale;
    p->getForce()[2] += Vp * mu_0 * modGradH[2] * _scale;
  }
}
*/

}
#endif
