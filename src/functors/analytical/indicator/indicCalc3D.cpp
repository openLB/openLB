/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Albert Mink
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

#include "indicCalc3D.h"
#include "indicCalc3D.hh"


namespace olb {

// arithmetic helper class for indicator 3d functors
template class IndicCalc3D<double,util::plus>;
template class IndicCalc3D<double,util::minus>;
template class IndicCalc3D<double,util::multiplies>;

template std::shared_ptr<IndicatorF3D<double>> operator+(std::shared_ptr<IndicatorF3D<double>> lhs,
    std::shared_ptr<IndicatorF3D<double>> rhs);
template std::shared_ptr<IndicatorF3D<double>> operator-(std::shared_ptr<IndicatorF3D<double>> lhs,
    std::shared_ptr<IndicatorF3D<double>> rhs);
template std::shared_ptr<IndicatorF3D<double>> operator*(std::shared_ptr<IndicatorF3D<double>> lhs,
    std::shared_ptr<IndicatorF3D<double>> rhs);

template std::shared_ptr<IndicatorF3D<double>> operator+(IndicatorIdentity3D<double> & lhs,
    std::shared_ptr<IndicatorF3D<double>>);
template std::shared_ptr<IndicatorF3D<double>> operator-(IndicatorIdentity3D<double> & lhs,
    std::shared_ptr<IndicatorF3D<double>>);
template std::shared_ptr<IndicatorF3D<double>> operator*(IndicatorIdentity3D<double> & lhs,
    std::shared_ptr<IndicatorF3D<double>>);

}

