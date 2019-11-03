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

/** \file
 * Interface for post-processing steps -- header file.
 */
#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include <vector>
#include "spatiallyExtendedObject2D.h"
#include "spatiallyExtendedObject3D.h"
#include "io/ostreamManager.h"

namespace olb {

/////////////////// Forward Declarations /////////////////////////////

template<typename T, typename DESCRIPTOR>
class BlockLatticeStructure2D;

template<typename T, typename DESCRIPTOR>
class BlockLatticeStructure3D;

template<typename T, typename DESCRIPTOR>
class BlockLattice2D;

template<typename T, typename DESCRIPTOR>
class BlockLattice3D;


/////////////////// 2D Postprocessing ///////////////////////////////

/// Interface of 2D post-processing steps.
template<typename T, typename DESCRIPTOR>
struct PostProcessor2D {
  virtual ~PostProcessor2D() { }
  /// Execute post-processing step
  virtual void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) =0;
  /// Execute post-processing step on a sublattice
  virtual void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_) =0;
  /// Extent of application area (0 for purely local operations)
  virtual int extent() const =0;
  /// Extent of application area along a direction (0 or 1)
  virtual int extent(int direction) const =0;
};

template<typename T, typename DESCRIPTOR>
class PostProcessorGenerator2D {
public:
  PostProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_);
  virtual ~PostProcessorGenerator2D() { }
  void shift(int deltaX, int deltaY);
  bool extract(int x0_, int x1_, int y0_, int y1_);
  virtual PostProcessor2D<T,DESCRIPTOR>* generate() const =0;
  virtual PostProcessorGenerator2D<T,DESCRIPTOR>* clone() const =0;
protected:
  int x0, x1, y0, y1;
};

template<typename T, typename DESCRIPTOR>
class LatticeCouplingGenerator2D {
public:
  LatticeCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_);
  virtual ~LatticeCouplingGenerator2D() { }
  void shift(int deltaX, int deltaY);
  bool extract(int x0_, int x1_, int y0_, int y1_);
  void reset(int x0_, int x1_, int y0_, int y1_);
  virtual PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D*> partners) const =0;
  virtual LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const =0;
protected:
  int x0, x1, y0, y1;
};


template<typename T, typename DESCRIPTOR>
struct LocalPostProcessor2D : public PostProcessor2D<T,DESCRIPTOR> {
};

template<typename T, typename DESCRIPTOR>
struct GlobalPostProcessor2D : public PostProcessor2D<T,DESCRIPTOR> {
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override =0;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_ ) override
  {
    this -> process(blockLattice);
  }
  int extent() const override
  {
    return 0;
  }
  int extent(int direction) const override
  {
    return 0;
  }
};


/////////////////// 3D Postprocessing ///////////////////////////////

template<typename T, typename DESCRIPTOR>
struct PostProcessor3D {
  virtual ~PostProcessor3D() { }
  /// Execute post-processing step
  virtual void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) =0;
  /// Execute post-processing step on a sublattice
  virtual void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_,
                                int z0_, int z1_ ) =0;
  /// Extent of application area (0 for purely local operations)
  virtual int extent() const =0;
  /// Extent of application area along a direction (0 or 1)
  virtual int extent(int direction) const =0;
};

template<typename T, typename DESCRIPTOR>
class PostProcessorGenerator3D {
public:
  PostProcessorGenerator3D( int x0_, int x1_, int y0_, int y1_,
                            int z0_, int z1_ );
  virtual ~PostProcessorGenerator3D() { }
  void shift(int deltaX, int deltaY, int deltaZ);
  bool extract(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  virtual PostProcessor3D<T,DESCRIPTOR>* generate() const =0;
  virtual PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const =0;
protected:
  int x0, x1, y0, y1, z0, z1;
};


template<typename T, typename DESCRIPTOR>
class LatticeCouplingGenerator3D {
public:
  LatticeCouplingGenerator3D( int x0_, int x1_, int y0_, int y1_,
                              int z0_, int z1_ );
  virtual ~LatticeCouplingGenerator3D() { }
  void shift(int deltaX, int deltaY, int deltaZ, int iC_=-1);
  bool extract(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  void reset(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  virtual PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject3D*> partners) const =0;
  virtual LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const =0;
protected:
  int x0, x1, y0, y1, z0, z1, iC;
};


template<typename T, typename DESCRIPTOR>
struct LocalPostProcessor3D : public PostProcessor3D<T,DESCRIPTOR> {
};

template<typename T, typename DESCRIPTOR>
struct GlobalPostProcessor3D : public PostProcessor3D<T,DESCRIPTOR> {
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override =0;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_,
                                int z0_, int z1_ ) override
  {
    this -> process(blockLattice);
  }
  int extent() const override
  {
    return 0;
  }
  int extent(int direction) const override
  {
    return 0;
  }
};

template<typename T, typename DESCRIPTOR>
struct StatisticsPostProcessor2D : public GlobalPostProcessor2D<T,DESCRIPTOR> {
  StatisticsPostProcessor2D();
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
};

template<typename T, typename DESCRIPTOR>
class StatPPGenerator2D
  : public PostProcessorGenerator2D<T,DESCRIPTOR> {
public:
  StatPPGenerator2D();
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>* clone() const override;
};

template<typename T, typename DESCRIPTOR>
struct StatisticsPostProcessor3D : public GlobalPostProcessor3D<T,DESCRIPTOR> {
  StatisticsPostProcessor3D();
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
};

template<typename T, typename DESCRIPTOR>
class StatPPGenerator3D
  : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  StatPPGenerator3D();
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const override;
};


}  // namespace olb

#endif
