/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Lukas Baron, Mathias J. Krause
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

#include "io/stlReader.h"
#include "indicatorF2D.h"
#include "indicatorF2D.hh"

namespace olb {

// stl indicator functors
template class STLreader<double>;

// indicator functors 2D
//template class IndicatorF2DfromIndicatorF3D<double>;
template class IndicatorCuboid2D<double>;
template class IndicatorCircle2D<double>;
template class IndicatorLayer2D<double>;
template IndicatorCuboid2D<double>* createIndicatorCuboid2D(XMLreader const& params, bool verbose);
template IndicatorCircle2D<double>* createIndicatorCircle2D(XMLreader const& params, bool verbose);

}
