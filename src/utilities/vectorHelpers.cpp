/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Lukas Baron, Mathias J. Krause
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

#include "vectorHelpers.h"
#include<vector>
#include<string>
#include<sstream>

namespace olb {
namespace util {

inline bool nearZero(const float& a);
inline bool nearZero(const double& a);


std::vector<double> fromVector3(const Vector<double, 3>& vec);
std::vector<double> fromVector2(const Vector<double, 2>& vec);

//void crossProduct3D(double c[3], const double a[3], const double b[3]);

inline void copyN(int c[], const int a[], const unsigned& dim);
inline void copyN(double c[], const double a[], const unsigned& dim);

inline void copy3(int c[], const int a[]);
inline void copy3(double c[], const double a[]);

std::vector<double> norm(const std::vector<double>& a);
std::vector<double> normalize(const std::vector<double>& a);

//bool triangleIntersection(const std::vector<double>& point0, const std::vector<double>& point1, const std::vector<double>& point2, const std::vector<double>& origign, const std::vector<double>& direction, double& distance);

//bool triangleIntersectionWithNormalDirection(const std::vector<double>& point0, const std::vector<double>& point1, const std::vector<double>& point2, const std::vector<double>& origign, const std::vector<double>& normalDirection, double& distance);


}
}
