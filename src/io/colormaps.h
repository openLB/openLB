/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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

#ifndef COLORMAPS_H
#define COLORMAPS_H

#include <string>
#include <vector>

namespace olb {

namespace graphics {

template <typename T>
struct ScalarFunction {
  virtual ~ScalarFunction() { }
  virtual T operator() (T x) const =0;
  virtual ScalarFunction<T>* clone() const=0;
};

template <typename T>
class LinearFunction : public ScalarFunction<T> {
public:
  LinearFunction(T x1_, T x2_, T y1_, T y2_);
  T operator() (T x) const override;
  LinearFunction<T>* clone() const override;
private:
  T x1, x2, y1, y2;
};

template <typename T>
class PowerLawFunction : public ScalarFunction<T> {
public:
  PowerLawFunction(T x1_, T x2_, T y1_, T y2_, T b_);
  T operator() (T x) const override;
  PowerLawFunction<T>* clone() const override;
private:
  T x1, x2, y1, y2;
  T b;
};

template <typename T>
struct Piece {
  Piece(T closedBegin_, T openEnd_) : closedBegin(closedBegin_), openEnd(openEnd_)
  {
  }
  T closedBegin, openEnd;
};

template <typename T>
class PiecewiseFunction : public ScalarFunction<T> {
public:
  PiecewiseFunction() { }
  ~PiecewiseFunction() override;
  PiecewiseFunction(PiecewiseFunction<T> const& rhs);
  PiecewiseFunction<T>& operator=(PiecewiseFunction<T> const& rhs);
  void swap(PiecewiseFunction<T>& rhs);
  void addPiece(Piece<T> piece, ScalarFunction<T>* f);
  T operator() (T x) const override;
  PiecewiseFunction<T>* clone() const override;
private:
  std::vector<Piece<T> > pieces;
  std::vector<ScalarFunction<T>*> functions;
};

template <typename T>
struct rgb {
  rgb(T r_, T g_, T b_) : r(r_), g(g_), b(b_)
  {
  }
  T r,g,b;
};

template <typename T>
class ColorMap {
public:
  ColorMap(PiecewiseFunction<T> const& red_,
           PiecewiseFunction<T> const& green_,
           PiecewiseFunction<T> const& blue_);
  rgb<T> get(T x) const;
private:
  PiecewiseFunction<T> red, green, blue;
};

namespace mapGenerators {

template <typename T> PiecewiseFunction<T> generateEarthRed();
template <typename T> PiecewiseFunction<T> generateEarthGreen();
template <typename T> PiecewiseFunction<T> generateEarthBlue();

template <typename T> PiecewiseFunction<T> generateWaterRed();
template <typename T> PiecewiseFunction<T> generateWaterGreen();
template <typename T> PiecewiseFunction<T> generateWaterBlue();

template <typename T> PiecewiseFunction<T> generateAirRed();
template <typename T> PiecewiseFunction<T> generateAirGreen();
template <typename T> PiecewiseFunction<T> generateAirBlue();

template <typename T> PiecewiseFunction<T> generateFireRed();
template <typename T> PiecewiseFunction<T> generateFireGreen();
template <typename T> PiecewiseFunction<T> generateFireBlue();

template <typename T> PiecewiseFunction<T> generateLeeLooRed();
template <typename T> PiecewiseFunction<T> generateLeeLooGreen();
template <typename T> PiecewiseFunction<T> generateLeeLooBlue();

template <typename T> ColorMap<T> generateMap(std::string mapName);

} // namespace mapGenerators

}  // namespace graphics

}  // namespace olb

#endif
