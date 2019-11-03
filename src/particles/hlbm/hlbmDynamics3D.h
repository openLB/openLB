/*  DESCRIPTOR Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2016 Thomas Henn, Fabian Klemens, Robn Trunk, Davide Dapelo
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

#ifndef PARTICLEDYNAMICS_3D_H
#define PARTICLEDYNAMICS_3D_H

#include "core/superLattice3D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF3D.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
class ParticleDynamics3D {
private:
  SuperLattice3D<T, DESCRIPTOR>& _sLattice;
  UnitConverter<T,DESCRIPTOR> const& _converter;
  SuperGeometry3D<T>& _superGeometry;
  std::shared_ptr<IndicatorF3D<T> > _indicatorF;
  std::vector<SmoothIndicatorF3D<T,T,true>* > _vectorOfIndicator;
  T _lengthX;
  T _lengthY;
  T _lengthZ;
  Vector<T,3> _accExt;
  bool _escapeFromDomain=false;
  bool _oldWallCollision=true;

public:
  ParticleDynamics3D(SuperLattice3D<T, DESCRIPTOR>& sLattice, 
                    UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperGeometry3D<T>& superGeometry, 
                    T lengthX, T lengthY, T lengthZ, 
                    Vector<T,3> accExt = Vector<T,3> (0.,0.,0.))
    : _sLattice(sLattice), 
      _converter(converter),
      _superGeometry(superGeometry), 
      _lengthX(lengthX), 
      _lengthY(lengthY), 
      _lengthZ(lengthZ),
      _accExt(accExt) 
  {}

  ParticleDynamics3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
		    UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperGeometry3D<T>& superGeometry, 
                    std::shared_ptr<IndicatorF3D<T> > indicatorF,
		    bool escapeFromDomain, bool oldWallCollision=false,
                    Vector<T,3> accExt = Vector<T,3> (0.,0.,0.))
      : _sLattice(sLattice),
	_converter(converter),
        _superGeometry(superGeometry),
	_indicatorF(indicatorF),
        _accExt(accExt),
	_escapeFromDomain(escapeFromDomain),
        _oldWallCollision(oldWallCollision)
  {}

  void addCuboid(Vector< T, 3> center, T xLength, T yLength, T zLength, T mass, T epsilon, Vector< T, 3 > theta, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  void addSphere(Vector< T, 3> center, T radius, T epsilon, T density, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  void addParticle(SmoothIndicatorF3D<T, T, true>& indicator);

  void computeBoundaryForce(std::vector<SmoothIndicatorF3D<T, T, true>* >& indicator);

  void addWallColl(SmoothIndicatorF3D<T, T, true>& indicator, T delta);

  void verletIntegration(SmoothIndicatorF3D<T, T, true>& indicator);

  void updateParticleDynamics(std::string name, SmoothIndicatorF3D<T, T, true>& indicator);

  void checkAndRemoveEscaped();

  void addParticleField(SmoothIndicatorF3D<T, T, true>& indicator);

  void simulateTimestep(std::string name);

  void print();

  void load(std::string filename, T epsilon);

  void save(std::string filename);

  void eulerIntegration(SmoothIndicatorF3D<T,T,true>& indicator);

};

}

#endif
