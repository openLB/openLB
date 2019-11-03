/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2015 Mathias J. Krause
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
 * Groups all the 2D include files in the complexDynamics directory.
*/

#include "advectionDiffusionDynamics.h"
#include "advectionDiffusionMomenta.h"
#include "ADMlatticeDescriptors.h"
#include "chopardDynamics.h"
#include "dynamics.h"
#include "dynOmegaLatticeDescriptors.h"
#include "dynSmagorinskyLatticeDescriptors.h"
#include "entropicDynamics.h"
#include "entropicLbHelpers.h"
#include "freeEnergyDynamics.h"
#include "freeEnergyPostProcessor2D.h"
#include "guoZhaoLatticeDescriptors.h"
#include "guoZhaoLbHelpers.h"
#include "guoZhaoDynamics.h"
#include "interactionPotential.h"
#include "latticeDescriptors.h"
#include "mrtDynamics.h"
#include "mrtLatticeDescriptors.h"
#include "navierStokesAdvectionDiffusionCouplingPostProcessor2D.h"
#include "navierStokesAdvectionDiffusionMRTCouplingPostProcessor2D.h"
#include "porousBGKdynamics.h"
#include "porousLatticeDescriptors.h"
#include "powerLawBGKdynamics.h"
#include "rtlbmDescriptors.h"
#include "shanChenDynOmegaForcedPostProcessor2D.h"
#include "shanChenForcedLatticeDescriptors.h"
#include "shanChenForcedPostProcessor2D.h"
#include "shanChenForcedSingleComponentPostProcessor2D.h"
#include "shearSmagorinskyLatticeDescriptors.h"
#include "smagorinskyBGKdynamics.h"
#include "smagorinskyGuoZhaoDynamics.h"
#include "SmagorinskyPorousParticleBGKdynamics.h"
#include "SmagorinskyPowerLawPorousBGKdynamics.h"
#include "SmagorinskyPowerLawBGKdynamics.h"
#include "smagorinskyMRTdynamics.h"
#include "stochasticSGSdynamics.h"
#include "superGuoZhaoPostProcessor2D.h"
