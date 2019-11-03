/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2015 Gilles Zahnd, Mathias J. Krause,
 *  Lukas Baron, Marie-Luise Maier
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

#include "frameChangeF3D.h"
#include "frameChangeF3D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class RotatingLinear3D<double>;
template class RotatingQuadratic1D<double>;
template class CirclePowerLaw3D<double>;
template class CirclePowerLawTurbulent3D<double>;
template class CirclePoiseuille3D<double>;
template class CirclePoiseuilleStrainRate3D<double, descriptors::D3Q19<>>;
template class RectanglePoiseuille3D<double>;
template class EllipticPoiseuille3D<double>;

template class AnalyticalPorousVelocity3D<double>;

template class AngleBetweenVectors3D<double,double>;

template class RotationRoundAxis3D<double,double>;

template class CylinderToCartesian3D<double,double>;

template class CartesianToCylinder3D<double,double>;

template class SphericalToCartesian3D<double,double>;

template class CartesianToSpherical3D<double,double>;

template class MagneticForceFromCylinder3D<double,double>;

template class MagneticFieldFromCylinder3D<double,double>;

}
