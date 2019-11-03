/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Marie-Luise Maier, Mathias J. Krause
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

/// Magnetic field that creates magnetization in wire has to be orthogonal to the wire.
/// to calculate the magnetic force on particles around a cylinder
/// (J. Lindner et al. / Computers and Chemical Engineering 54 (2013) 111-121)

#include "dynamics/latticeDescriptors.h"
 
#include "particles/particle3D.h"
#include "particles/particle3D.hh"
#include "interpMagForceForMagP3D.h"
#include "interpMagForceForMagP3D.hh"

namespace olb {

template class InterpMagForceForMagP3D<double, MagneticParticle3D, descriptors::D3Q19<>>;

}  // namespace olb

