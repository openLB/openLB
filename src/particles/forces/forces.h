/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause
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
 * Groups all the include .h-files for 3D particles forces in
 * the particles forces directory.
 */

#include "force3D.h"
#include "stokesDragForce3D.h"
#include "weightForce3D.h"
#include "buoyancyForce3D.h"
#include "hertzMindlinDeresiewicz3D.h"
#include "transferExternalForce3D.h"
#include "interpMagForceForMagP3D.h"
#include "magneticForceForMagP3D.h"
#include "linearDampingForceForMagDipoleMoment3D.h"
#include "stokesDragForceForHomVelField3D.h"
#include "magneticForceFromHField3D.h"
#include "forceFromExtField3D.h"


