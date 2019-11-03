#  This file is part of the OpenLB library
#
#  Copyright (C) 2007, 2017 Markus Mohrhard, Mathias Krause
#  E-mail contact: info@openlb.net
#  The most recent release of OpenLB can be downloaded at
#  <http://www.openlb.net/>
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public 
#  License along with this program; if not, write to the Free 
#  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#  Boston, MA  02110-1301, USA.

###########################################################################
###########################################################################

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))

include $(dir $(mkfile_path))/config.mk

###########################################################################
## conditional settings

ifeq ($(BUILDTYPE), precompiled)
   CXXFLAGS := -DOLB_PRECOMPILED $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), MPI)
   CXXFLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS) $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), OMP)
   CXXFLAGS := -DPARALLEL_MODE_OMP $(OMPFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

ifeq ($(PARALLEL_MODE), HYBRID)
   CXXFLAGS := -DPARALLEL_MODE_OMP -DPARALLEL_MODE_MPI $(OMPFLAGS) $(MPIFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

LDFLAGS += $(if $(filter $(FEATURES), OPENBLAS),-lopenblas)

CXXFLAGS += $(foreach feature,$(FEATURES),-DFEATURE_$(feature))

###########################################################################
## defines shell

SHELL           := /bin/sh

###########################################################################
## dependencies, object, library directory and library name

DEPENDDIR       := build/$(BUILDTYPE)/dep
OBJDIR          := build/$(BUILDTYPE)/obj
LIBDIR          := build/$(BUILDTYPE)/lib
LIB             := olb
LIBS            := -l$(LIB) -lz

###########################################################################
## search directories

SUBDIRS         := src/boundary \
                   src/communication \
                   src/dynamics \
                   src/core \
                   src/geometry \
                   src/external/tinyxml \
                   src/external/zlib \
                   src/functors \
                   src/functors/analytical \
                   src/functors/analytical/indicator \
                   src/functors/lattice \
                   src/functors/lattice/indicator \
                   src/functors/lattice/integral \
                   src/io \
                   src/particles \
                   src/particles/forces \
                   src/particles/boundaries \
                   src/utilities \

EXAMPLEDIRS     := examples/laminar/bstep2d \
                   examples/laminar/bstep3d \
                   examples/laminar/cavity2d/sequential \
                   examples/laminar/cavity2d/parallel \
                   examples/laminar/cavity3d/sequential \
                   examples/laminar/cavity3d/parallel \
                   examples/laminar/cylinder2d \
                   examples/laminar/cylinder3d \
                   examples/laminar/poiseuille2d \
                   examples/laminar/poiseuille3d \
                   examples/laminar/powerLaw2d \
                   examples/multiComponent/contactAngle2d \
                   examples/multiComponent/contactAngle3d \
                   examples/multiComponent/microFluidics2d \
                   examples/multiComponent/phaseSeparation2d \
                   examples/multiComponent/phaseSeparation3d \
                   examples/multiComponent/rayleighTaylor2d \
                   examples/multiComponent/rayleighTaylor3d \
                   examples/multiComponent/youngLaplace2d \
                   examples/multiComponent/youngLaplace3d \
                   examples/particles/bifurcation3d/eulerEuler \
                   examples/particles/bifurcation3d/eulerLagrange \
                   examples/particles/dkt2d \
                   examples/particles/magneticParticles3d \
                   examples/particles/settlingCube3d \
                   examples/porousMedia/porousPoiseuille2d \
                   examples/porousMedia/porousPoiseuille3d \
                   examples/thermal/porousPlate2d \
                   examples/thermal/porousPlate3d \
                   examples/thermal/rayleighBenard2d \
                   examples/thermal/rayleighBenard3d \
                   examples/thermal/squareCavity2d \
                   examples/thermal/squareCavity3d \
                   examples/turbulence/aorta3d \
                   examples/turbulence/nozzle3d \
                   examples/turbulence/tgv3d \
                   examples/turbulence/venturi3d \

INCLUDEDIRS     := src \
                   src/ \
                   src/external \
                   src/external/zlib

BUILDTYPEDIRS   := build/precompiled \
                   build/generic

IDIR            := $(foreach d,$(INCLUDEDIRS),-I$(ROOT)/$(d))
