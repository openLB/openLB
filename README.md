# OpenLB – Open Source Lattice Boltzmann Code
Unofficial mirror

The OpenLB project provides a C++ package for the implementation of lattice Boltzmann methods that is general enough to address a vast range of tansport problems, e.g. in computational fluid dynamics. The source code is publicly available and constructed in a well readable, modular way. This enables for a fast implementation of both simple academic test problems and advanced engineering applications. It is also easily extensible to include new physical content.

## New user-friendly features:
* New meta descriptor concept
* New homogenised lattice Boltzmann method
* New free energy model
* Validated wall shear stress functor

## New examples:
* multiComponent/contactAngle2d and multiComponent/contactAngle3d
* multiComponent/youngLaplace2d and multiComponent/youngLaplace3d
* multiComponent/microFluidics2d
* particles/magneticParticles3d
* particles/settlingCube3d
* particles/dkt2d
* porousMedia/porousPoiseuille3d

## Minor improvements and developer notes:
* Restructure example folder
* Restructure and improve functors
* New std::shared_ptr-based functor arithmetic to ease memory management and enable functor composition
* Convenient relative and absolute Lp error norm functors
* Bug fixed in GnuPlotWriter
* C++ 14 standard is now mandatory
