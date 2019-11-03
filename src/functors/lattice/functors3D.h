/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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
 * Groups all include files for the directory genericFunctions.
 */

#include "blockBaseF2D.h"
#include "blockBaseF3D.h"
#include "blockCalcF3D.h"
#include "blockLatticeIntegralF3D.h"
#include "blockMin3D.h"
#include "blockMax3D.h"
#include "blockAverage3D.h"
#include "blockGeometryFaces3D.h"
#include "blockLatticeLocalF3D.h"
#include "blockLocalAverage3D.h"
#include "latticeFrameChangeF3D.h"
#include "reductionF3D.h"
#include "blockReduction3D2D.h"
#include "superBaseF3D.h"
#include "superCalcF3D.h"
#include "superConst3D.h"
#include "superLatticeIntegralF3D.h"
#include "superMin3D.h"
#include "superMax3D.h"
#include "superAverage3D.h"
#include "superGeometryFaces3D.h"
#include "superLatticeLocalF3D.h"
#include "superLocalAverage3D.h"
#include "turbulentF3D.h"
#include "superErrorNorm3D.h"
#include "indicator/indicator3D.h"
#include "integral/integral3D.h"
#include "timeAveraged/superLatticeTimeAveraged2D.h"
