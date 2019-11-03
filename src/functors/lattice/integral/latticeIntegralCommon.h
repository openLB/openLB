/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef LATTICE_INTEGRAL_COMMON_H
#define LATTICE_INTEGRAL_COMMON_H

namespace olb {


/// Lp norm functor implementation details specific to the P parameter
/**
 * Used in BlockLpNorm2D, BlockLpNorm3D, SuperLpNorm2D and SuperLpNorm3D.
 **/
template <typename T, typename W, int P>
struct LpNormImpl {
  inline W operator()(W output, W tmp, T weight);
  inline W enclose(W output);
};

template <typename T, typename W, int P>
inline W LpNormImpl<T,W,P>::operator()(W output, W tmp, T weight)
{
  return output + pow(fabs(tmp), P)*weight;
}

template <typename T, typename W, int P>
inline W LpNormImpl<T,W,P>::enclose(W output)
{
  return pow(output, 1. / P);
}

/// Linf norm functor implementation details
template <typename T, typename W>
struct LpNormImpl<T,W,0> {
  inline W operator()(W output, W tmp, T weight)
  {
    return std::max(output, fabs(tmp));
  }
  inline W enclose(W output)
  {
    return output;
  }
};

/// L1 norm functor implementation details
template <typename T, typename W>
struct LpNormImpl<T,W,1> {
  inline W operator()(W output, W tmp, T weight)
  {
    return output + fabs(tmp)*weight;
  }
  inline W enclose(W output)
  {
    return output;
  }
};

/// L2 norm functor implementation details
template <typename T, typename W>
struct LpNormImpl<T,W,2> {
  inline W operator()(W output, W tmp, T weight)
  {
    return output + tmp*tmp*weight;
  }
  inline W enclose(W output)
  {
    return sqrt(output);
  }
};


} // end namespace olb

#endif
