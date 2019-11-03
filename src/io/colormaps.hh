/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#ifndef COLORMAPS_HH
#define COLORMAPS_HH

#include "colormaps.h"
#include <cmath>

namespace olb {

namespace graphics {

template <typename T>
LinearFunction<T>::LinearFunction(T x1_, T x2_, T y1_, T y2_)
  : x1(x1_), x2(x2_), y1(y1_), y2(y2_)
{ }

template <typename T>
T LinearFunction<T>::operator() (T x) const
{
  return ( (y2-y1) * x + x2*y1-x1*y2 )/(x2-x1);
}

template <typename T>
LinearFunction<T>* LinearFunction<T>::clone() const
{
  return new LinearFunction(*this);
}

template <typename T>
PowerLawFunction<T>::PowerLawFunction(T x1_, T x2_, T y1_, T y2_, T b_)
  : x1(std::pow(x1_,b_)), x2(std::pow(x2_,b_)), y1(y1_), y2(y2_), b(b_)
{ }

template <typename T>
T PowerLawFunction<T>::operator() (T x) const
{
  return ( (y2-y1) * std::pow(x,b) + x2*y1-x1*y2 )/(x2-x1);
}

template <typename T>
PowerLawFunction<T>* PowerLawFunction<T>::clone() const
{
  return new PowerLawFunction(*this);
}


template <typename T>
PiecewiseFunction<T>::PiecewiseFunction(PiecewiseFunction<T> const& rhs)
  : pieces(rhs.pieces),
    functions(rhs.functions.size())
{
  for (unsigned iF=0; iF<functions.size(); ++iF) {
    functions[iF] = rhs.functions[iF]->clone();
  }
}

template <typename T>
PiecewiseFunction<T>& PiecewiseFunction<T>::operator=(PiecewiseFunction<T> const& rhs)
{
  PiecewiseFunction<T>(rhs).swap(*this);
  return *this;
}

template <typename T>
PiecewiseFunction<T>::~PiecewiseFunction()
{
  for (unsigned iF=0; iF<functions.size(); ++iF) {
    delete functions[iF];
  }
}

template <typename T>
void PiecewiseFunction<T>::swap(PiecewiseFunction<T>& rhs)
{
  pieces.swap(rhs.pieces);
  functions.swap(rhs.functions);
}

template <typename T>
void PiecewiseFunction<T>::addPiece(Piece<T> piece, ScalarFunction<T>* f)
{
  typename std::vector<Piece<T> >::iterator pieceIt = pieces.begin();
  typename std::vector<ScalarFunction<T>*>::iterator fIt = functions.begin();
  while (pieceIt != pieces.end() && piece.closedBegin >= pieceIt->closedBegin) {
    ++pieceIt;
    ++fIt;
  }
  pieces.insert(pieceIt, piece);
  functions.insert(fIt, f);
}

template <typename T>
T PiecewiseFunction<T>::operator() (T x) const
{
  if (pieces.empty() || x<pieces[0].closedBegin) {
    return T();
  }
  unsigned iPiece=0;
  while (iPiece != pieces.size() && x >= pieces[iPiece].openEnd) {
    ++iPiece;
  }
  if (iPiece == pieces.size() || x < pieces[iPiece].closedBegin) {
    return T();
  }
  return (*functions[iPiece])(x);
}

template <typename T>
PiecewiseFunction<T>* PiecewiseFunction<T>::clone() const
{
  return new PiecewiseFunction(*this);
}

template <typename T>
ColorMap<T>::ColorMap(PiecewiseFunction<T> const& red_,
                      PiecewiseFunction<T> const& green_,
                      PiecewiseFunction<T> const& blue_)
  : red(red_), green(green_), blue(blue_)
{ }

template <typename T>
rgb<T> ColorMap<T>::get(T x) const
{
  return rgb<T>( red(x), green(x), blue(x) );
}

namespace mapGenerators {

template <typename T>
PiecewiseFunction<T> generateEarthRed()
{
  T p0  = (T) 0.;
  T p1  = (T) 3./8.;
  T p2  = (T) 6./8.;
  T p3  = (T) 1.;

  PiecewiseFunction<T> earthRed;
  earthRed.addPiece (Piece<T>(p0, p1), new PowerLawFunction<T>(p0, p1, (T)0., (T)0.8, (T)0.6) );
  earthRed.addPiece (Piece<T>(p1, p2), new PowerLawFunction<T>(p1, p2, (T)0.8, (T)0.9, (T)0.9) );
  earthRed.addPiece (Piece<T>(p2, p3), new PowerLawFunction<T>(p2, p3, (T)0.9, (T)1.0, (T)0.2) );

  return earthRed;
}

template <typename T>
PiecewiseFunction<T> generateEarthGreen()
{
  T p0  = (T) 0.;
  T p1  = (T) 3./8.;
  T p2  = (T) 6./8.;
  T p3  = (T) 1.;

  PiecewiseFunction<T> earthGreen;
  earthGreen.addPiece (Piece<T>(p0, p1), new PowerLawFunction<T>(p0, p1, (T)0., (T)0.5, (T)0.6) );
  earthGreen.addPiece (Piece<T>(p1, p2), new PowerLawFunction<T>(p1, p2, (T)0.5, (T)0.9, (T)0.2) );
  earthGreen.addPiece (Piece<T>(p2, p3), new PowerLawFunction<T>(p2, p3, (T)0.9, (T)1.0, (T)0.2) );

  return earthGreen;
}

template <typename T>
PiecewiseFunction<T> generateEarthBlue()
{
  T p0  = (T) 0.;
  T p1  = (T) 3./8.;
  T p2  = (T) 6./8.;
  T p3  = (T) 1.;

  PiecewiseFunction<T> earthBlue;
  earthBlue.addPiece (Piece<T>(p0, p1), new PowerLawFunction<T>(p0, p1, (T)0., (T)0.5, (T)0.6) );
  earthBlue.addPiece (Piece<T>(p1, p2), new PowerLawFunction<T>(p1, p2, (T)0.5, (T)0.7, (T)0.2) );
  earthBlue.addPiece (Piece<T>(p2, p3), new PowerLawFunction<T>(p2, p3, (T)0.7, (T)1.0, (T)0.2) );

  return earthBlue;
}

template <typename T>
PiecewiseFunction<T> generateWaterRed()
{
  T p0  = (T) 0.;
  T p1  = (T) 3./8.;
  T p2  = (T) 6./8.;
  T p3  = (T) 1.;

  PiecewiseFunction<T> waterRed;
  waterRed.addPiece (Piece<T>(p0, p1), new PowerLawFunction<T>(p0, p1, (T)0., (T)0.5, (T)0.6) );
  waterRed.addPiece (Piece<T>(p1, p2), new PowerLawFunction<T>(p1, p2, (T)0.5, (T)0.7, (T)0.2) );
  waterRed.addPiece (Piece<T>(p2, p3), new PowerLawFunction<T>(p2, p3, (T)0.7, (T)1.0, (T)0.2) );

  return waterRed;
}

template <typename T>
PiecewiseFunction<T> generateWaterGreen()
{
  T p0  = (T) 0.;
  T p1  = (T) 3./8.;
  T p2  = (T) 6./8.;
  T p3  = (T) 1.;

  PiecewiseFunction<T> waterGreen;
  waterGreen.addPiece (Piece<T>(p0, p1), new PowerLawFunction<T>(p0, p1, (T)0., (T)0.5, (T)0.6) );
  waterGreen.addPiece (Piece<T>(p1, p2), new PowerLawFunction<T>(p1, p2, (T)0.5, (T)0.9, (T)0.2) );
  waterGreen.addPiece (Piece<T>(p2, p3), new PowerLawFunction<T>(p2, p3, (T)0.9, (T)1.0, (T)0.2) );

  return waterGreen;
}

template <typename T>
PiecewiseFunction<T> generateWaterBlue()
{
  T p0  = (T) 0.;
  T p1  = (T) 3./8.;
  T p2  = (T) 6./8.;
  T p3  = (T) 1.;

  PiecewiseFunction<T> waterBlue;
  waterBlue.addPiece (Piece<T>(p0, p1), new PowerLawFunction<T>(p0, p1, (T)0., (T)0.8, (T)0.6) );
  waterBlue.addPiece (Piece<T>(p1, p2), new PowerLawFunction<T>(p1, p2, (T)0.8, (T)0.9, (T)0.9) );
  waterBlue.addPiece (Piece<T>(p2, p3), new PowerLawFunction<T>(p2, p3, (T)0.9, (T)1.0, (T)0.2) );

  return waterBlue;
}

template <typename T>
PiecewiseFunction<T> generateAirRed()
{
  T p0  = (T) 0.;
  T p1  = (T) 1.;

  PiecewiseFunction<T> airRed;
  airRed.addPiece (Piece<T>(p0, p1), new LinearFunction<T>(p0, p1, (T)0., (T)1.) );

  return airRed;
}

template <typename T>
PiecewiseFunction<T> generateAirGreen()
{
  T p0  = (T) 0.;
  T p1  = (T) 1.;

  PiecewiseFunction<T> airGreen;
  airGreen.addPiece (Piece<T>(p0, p1), new LinearFunction<T>(p0, p1, (T)1., (T)0.) );

  return airGreen;
}

template <typename T>
PiecewiseFunction<T> generateAirBlue()
{
  T p0  = (T) 0.;
  T p1  = (T) 1.;

  PiecewiseFunction<T> airBlue;
  airBlue.addPiece (Piece<T>(p0, p1), new LinearFunction<T>(p0, p1, (T)1., (T)1.) );

  return airBlue;
}


template <typename T>
PiecewiseFunction<T> generateFireRed()
{
  T p0  = (T) 0.;
  T p1  = (T) 0.36;
  T p3  = (T) 1.;

  PiecewiseFunction<T> fireRed;
  fireRed.addPiece (Piece<T>(p0, p1), new LinearFunction<T>(p0, p1, (T)0., (T)1.) );
  fireRed.addPiece (Piece<T>(p1, p3), new LinearFunction<T>(p1, p3, (T)1., (T)1.) );

  return fireRed;
}

template <typename T>
PiecewiseFunction<T> generateFireGreen()
{
  T p0  = (T) 0.;
  T p1  = (T) 0.36;
  T p2  = (T) 0.75;
  T p3  = (T) 1.;

  PiecewiseFunction<T> fireGreen;
  fireGreen.addPiece (Piece<T>(p0, p1), new LinearFunction<T>(p0, p1, (T)0., (T)0.) );
  fireGreen.addPiece (Piece<T>(p1, p2), new LinearFunction<T>(p1, p2, (T)0., (T)1.) );
  fireGreen.addPiece (Piece<T>(p2, p3), new LinearFunction<T>(p2, p3, (T)1., (T)1.) );

  return fireGreen;
}

template <typename T>
PiecewiseFunction<T> generateFireBlue()
{
  T p0  = (T) 0.;
  T p2  = (T) 0.75;
  T p3  = (T) 1.;

  PiecewiseFunction<T> fireBlue;
  fireBlue.addPiece (Piece<T>(p0, p2), new LinearFunction<T>(p0, p2, (T)0., (T)0.) );
  fireBlue.addPiece (Piece<T>(p2, p3), new LinearFunction<T>(p2, p3, (T)0., (T)1.) );

  return fireBlue;
}


template <typename T>
PiecewiseFunction<T> generateLeeLooRed()
{
  T p0  = (T) 0.;
  T p2  = (T) 3./8.;
  T p3  = (T) 5./8.;
  T p4  = (T) 7./8.;
  T p5  = (T) 1.;
  T p6  = (T) 9./8.;


  PiecewiseFunction<T> leeLooRed;
  leeLooRed.addPiece (Piece<T>(p0, p2), new LinearFunction<T>(p0, p2, (T)0., (T)0.) );
  leeLooRed.addPiece (Piece<T>(p2, p3), new LinearFunction<T>(p2, p3, (T)0., (T)1.) );
  leeLooRed.addPiece (Piece<T>(p3, p4), new LinearFunction<T>(p3, p4, (T)1., (T)1.) );
  leeLooRed.addPiece (Piece<T>(p4, p5), new LinearFunction<T>(p4, p6, (T)1., (T)0.) );

  return leeLooRed;
}

template <typename T>
PiecewiseFunction<T> generateLeeLooGreen()
{
  T p0  = (T) 0.;
  T p1  = (T) 1./8.;
  T p2  = (T) 3./8.;
  T p3  = (T) 5./8.;
  T p4  = (T) 7./8.;
  T p5  = (T) 1.;
  T p6  = (T) 9/8;


  PiecewiseFunction<T> leeLooGreen;
  leeLooGreen.addPiece (Piece<T>(p0, p1), new LinearFunction<T>(p0, p1, (T)0., (T)0.) );
  leeLooGreen.addPiece (Piece<T>(p1, p2), new LinearFunction<T>(p1, p2, (T)0., (T)1.) );
  leeLooGreen.addPiece (Piece<T>(p2, p3), new LinearFunction<T>(p2, p3, (T)1., (T)1.) );
  leeLooGreen.addPiece (Piece<T>(p3, p4), new LinearFunction<T>(p3, p4, (T)1., (T)0.) );
  leeLooGreen.addPiece (Piece<T>(p4, p5), new LinearFunction<T>(p4, p6, (T)0., (T)0.) );

  return leeLooGreen;
}

template <typename T>
PiecewiseFunction<T> generateLeeLooBlue()
{
  T pm1 = (T) -1./8.;
  T p0  = (T) 0.;
  T p1  = (T) 1./8.;
  T p2  = (T) 3./8.;
  T p3  = (T) 5./8.;
  T p5  = (T) 1.;


  PiecewiseFunction<T> leeLooBlue;
  leeLooBlue.addPiece (Piece<T>(p0, p1), new LinearFunction<T>(pm1, p1, (T)0., (T)1.) );
  leeLooBlue.addPiece (Piece<T>(p1, p2), new LinearFunction<T>(p1, p2, (T)1., (T)1.) );
  leeLooBlue.addPiece (Piece<T>(p2, p3), new LinearFunction<T>(p2, p3, (T)1., (T)0.) );
  leeLooBlue.addPiece (Piece<T>(p3, p5), new LinearFunction<T>(p3, p5, (T)0., (T)0.) );

  return leeLooBlue;
}

template <typename T>
ColorMap<T> generateMap(std::string mapName)
{
  if (mapName == "earth") {
    return ColorMap<T> (
             generateEarthRed<T>(),
             generateEarthGreen<T>(),
             generateEarthBlue<T>() );
  } else if (mapName == "water") {
    return ColorMap<T> (
             generateWaterRed<T>(),
             generateWaterGreen<T>(),
             generateWaterBlue<T>() );
  } else if (mapName == "air") {
    return ColorMap<T> (
             generateAirRed<T>(),
             generateAirGreen<T>(),
             generateAirBlue<T>() );
  } else if (mapName == "fire") {
    return ColorMap<T> (
             generateFireRed<T>(),
             generateFireGreen<T>(),
             generateFireBlue<T>() );
  } else if (mapName == "leeloo") {
    return ColorMap<T> (
             generateLeeLooRed<T>(),
             generateLeeLooGreen<T>(),
             generateLeeLooBlue<T>() );
  }
  return ColorMap<T> (
           generateLeeLooRed<T>(),
           generateLeeLooGreen<T>(),
           generateLeeLooBlue<T>() );
}

} // namespace mapGenerators

}  // namespace graphics

}  // namespace olb

#endif
