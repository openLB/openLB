/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause
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

#ifndef SUPERPARTICLESYSTEM_3D_H
#define SUPERPARTICLESYSTEM_3D_H

#define shadows

#include <set>
#include <vector>
#include <list>
#include <deque>
#include <memory>
#include "boundary/boundary3D.h"
#include "communication/loadBalancer.h"
#include "communication/mpiManager.h"
#include "core/superLattice3D.h"
#include "forces/force3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "functors/analytical/analyticalF.h"
#include "geometry/cuboidGeometry3D.h"
#include "geometry/superGeometry3D.h"
#include "particleSystem3D.h"
#include "superParticleSysVTUout.h"
#include "functors/lattice/superLatticeLocalF3D.h"
#include "twoWayCouplings/twoWayCouplings3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class SuperParticleSysVtuWriter;

template<typename T>
class SuperParticleSysVtuWriterMag;

template<typename T, template<typename U> class PARTICLETYPE>
class Force3D;

template<typename T, template<typename U> class PARTICLETYPE>
class Boundary3D;

template <typename T, typename DESCRIPTOR>
class SuperLatticeInterpVelocity3D;

/**
 *  The class superParticleSystem is the basis for particulate flows within OpenLB.
 *  Use one of the constructors to instantiate a superParticleSystem. This creates
 *  one particleSystem for each cuboid according to the structure found in
 *  cuboid Geometry.
 *  <UL>
 *  <LI>Add single or several particles using one of the addParticle() functions.</LI>
 *  <LI>Add forces acting on the particles using the addForce() function.</LI>
 *  <LI>Add particle boundaries using the addBoundary() function.</LI>
 *  <LI>Finally compute one timestep using the simulate() function.</LI>
 *  </UL>
 */

template<typename T, template<typename U> class PARTICLETYPE>
class SuperParticleSystem3D : public SuperStructure3D<T> {

public:
  time_t _stopSorting;

  /// Constructor for SuperParticleSystem
  SuperParticleSystem3D(CuboidGeometry3D<T>& cuboidGeometry,
                        LoadBalancer<T>& loadBalancer, SuperGeometry3D<T>&);
  SuperParticleSystem3D(SuperGeometry3D<T>&);

  /// Copy Constructor for SuperParticleSystem
  SuperParticleSystem3D(SuperParticleSystem3D<T, PARTICLETYPE>& spSys);
  SuperParticleSystem3D(SuperParticleSystem3D<T, PARTICLETYPE> const& spSys);

  /// Move Constructor for SuperParticleSystem
  SuperParticleSystem3D(SuperParticleSystem3D<T, PARTICLETYPE> && spSys);
  /// Destructor
  ~SuperParticleSystem3D() override {};

  /// Add a Particle to SuperParticleSystem
  void addParticle(PARTICLETYPE<T> &p);

  /// Add a number of identical particles randomly distributed in a given IndicatorF3D
  void addParticle(IndicatorF3D<T>& ind, T mas, T rad, int no = 1, std::vector<T> vel = {0., 0., 0.});
  void addParticle(IndicatorF3D<T>& ind, T mas, T rad, int no, int id,
                   std::vector<T> vel, std::vector<T> dMoment, std::vector<T> aVel,
                   std::vector<T> torque, T magnetisation, int sActivity)
  {
    addParticle(ind, mas, rad, no, vel);
  };

  /// Add a number of identical Pprticles randomly distributed in a given IndicatorF3D
  /// and in given Material Number
  void addParticle(IndicatorF3D<T>& ind, std::set<int>  material, T mas, T rad, int no = 1,
                   std::vector<T> vel = {0., 0., 0.});
  void addParticle(IndicatorF3D<T>& ind,  std::set<int>  material, T mas, T rad, int no, int id,
                   std::vector<T> vel, std::vector<T> dMoment, std::vector<T> aVel,
                   std::vector<T> torque, T magnetisation, int sActivity)
  {
    addParticle(ind, material, mas, rad, no, vel);
  };

  /// Add a number of identical particles randomly distributed in a given Material Number
  void addParticle(std::set<int>  material, int no, T mas, T rad, std::vector<T> vel = {0., 0., 0.});

  /// Add a number of identical particles equally distributed in a given Material Number
  void addParticleEquallyDistributed(IndicatorCuboid3D<T>& cuboid, T pMass,
                                     T pRad,/* number of particles on x, y, z axis*/
                                     int nox, int noy, int noz, std::vector<T> vel = {0., 0., 0.});

  void addParticleEquallyDistributed(IndicatorCuboid3D<T>& cuboid, int nox, int noy, int noz, PARTICLETYPE<T>& p);

  /// Generates particle at a circle shaped inlet, amount given by mass concentration in feedstream.
  /// It is taken care that the particles do not overlap during initialization, therefore they are
  /// saved in posDeq with size deqSize.
  /// The parameters particlesPerPhyTimeStep and inletVec are initialized for iT = 0 and given back
  /// by reference. They can be defined by an arbitrary value.
  template<typename DESCRIPTOR>
  void generateParticlesCircleInletMassConcentration(
    IndicatorCircle3D<T>& indicatorCircle, T particleMassConcentration, T charPhysVelocity,
    T conversionFactorTime, SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel,
    PARTICLETYPE<T>& p, std::set<int> material, int iT, T& particlesPerPhyTimeStep,
    std::vector<T>& inletVec, std::deque<std::vector<T>>& posDeq, int deqSize);

  /// Add particles form a File. Save using saveToFile(std::string name)
  void addParticlesFromFile(std::string name, T mass, T radius);
  /// Add a number of particles with a certain ID (TracerParticle) equally distributed in a given IndicatorF3D
  void addTracerParticle(IndicatorF3D<T>& ind, T idTP, T mas, T rad, int noTP = 1, std::vector<T> vel = {0., 0., 0.});

  /// Add a number of unidentical particles with normally distributed radius (Box-Muller Method) in a given IndicatorF3D
  /// with specific appearance probability
  void addParticleBoxMuller(IndicatorF3D<T>& ind, T partRho, T mu, T sigma, int no = 1, T appProb = 1., std::vector<T> vel
                            = {0., 0., 0.});

  /// Removes all particles from System
  void clearParticles();

  /// Generates particles with specific volume concentration conc equally
  /// and randomly distributed in given IndicatorCuboid maintaining a minimum
  /// distance between each other.
  /// Be aware that long calculation time can occur because of minDist check.
  void addParticleWithDistance(IndicatorCuboid3D<T>& ind,
                               T pMass, T pRad, std::vector<T> vel,
                               T conc, // volume concentration of particles, noP*vol_1p/volF = conc
                               T minDist, // minimum distance between each particle
                               bool checkDist // check whether minDist is choosen too large
                              );

  /// Integrate on Timestep dT, scale = true keeps the particle velocity in stable range
  void simulate(T dT, bool scale = false);
  // multiple collision models
  void simulate(T dT, std::set<int> sActivityOfFreeParticle, bool scale = false)
  {
    simulate(dT, scale);
  };
  /// Integrate on Timestep dT with two-way coupling, scale = true keeps the particle velocity in stable range
  void simulateWithTwoWayCoupling_Mathias ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps = 1, bool resetExternalField = true, bool scale = false );
  void simulateWithTwoWayCoupling_Davide ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps = 1, bool resetExternalField = true, bool scale = false );

  /// Gives random velocity to all particles
  void setParticlesVelRandom(T velFactor);
  /// Changes particle positions randomly
  void setParticlesPosRandom(T posFactor);
  void setParticlesPosRandom(T posFactorX, T posFactorY, T posFactorZ);
  /// Gives random dipolemoment orientation to all MagneticParticle3D
  void setMagneticParticlesdMomRandom() {};
  /// Gives specific attributes to all MagneticParticle3D
  void setMagneticParticles(std::vector<T> dMoment, std::vector<T> vel, std::vector<T> aVel, std::vector<T> torque,
                            T magnetisation) {};
  void setMagneticParticles(std::vector<T> dMoment, std::vector<T> vel, std::vector<T> aVel, std::vector<T> torque,
                            T magnetisation, int sActivity) {};

  /// Agglomerate detection functions: Todo: enable for paralle mode
  /// Initializes an empty agglomerate list in every particleSystem3D
  void prepareAgglomerates() {};
  /// Adds new generated particles to the list of non agglomerated Particles
  void initAggloParticles() {};
  /// Detects and manages particle agglomerates
  void findAgglomerates(int iT, int itVtkOutputMagParticles) {};

  /// Tests if particles with specific sActivity exist
  bool particleSActivityTest(int sActivity)
  {
    return 0;
  };

  /// Set overlap of ParticleSystems, overlap has to be in lattice units
  /// particle system _overlap+1 <= _superGeometry.getOverlap()
  void setOverlap(T);
  /// Get overlap of ParticleSystems
  T getOverlap();

  /// Save Particles to file. Add using addParticlesFromFile(std::string name, T mass, T radius);
  void saveToFile(std::string name);

  /// Get global number of particles
  int globalNumOfParticles();
  /// Get global number of shadow particles (particles hold in overlap)
  int globalNumOfShadowParticles();
  /// Get global number of active particles
  int globalNumOfActiveParticles();
  /// Get number of particles computed on this node
  int rankNumOfParticles();
  /// Get number of shadow particles computed on this node
  int rankNumOfShadowParticles();
  /// Get number of active particles computed on this node
  int rankNumOfActiveParticles();
  /// Get number of TracerParticles computed on this node
  int rankNumOfTracerParticles();
  /// Get number of TracerParticles computed on this node
  int globalNumOfTracerParticles();

  /// Get ParticleSystems
  std::vector<ParticleSystem3D<T, PARTICLETYPE>*> getParticleSystems();
  /// Get ParticleSystems
  std::vector<ParticleSystem3D<T, PARTICLETYPE>*>& getPSystems();
  /// Get a ParticleSystem
  ParticleSystem3D<T, PARTICLETYPE>& operator[](int i);
  /// Get number of linked Forces
  std::vector<int> numOfForces();

  /// Get number of ParticleSystems
  int numOfPSystems();

  /// Get number of particles in the vicinity of material number mat
  int countMaterial(int mat);

  /// Add a force to system
  void addForce(std::shared_ptr<Force3D<T, PARTICLETYPE> > f);
  /// Add a boundary to system
  void addBoundary(std::shared_ptr<Boundary3D<T, PARTICLETYPE> > b);

  /// Set particle velocity to fluid velocity (e.g. as inital condition
  template<typename DESCRIPTOR>
  void setVelToFluidVel(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>&);

  /// Set particle velocity to analytical velocity (e.g. as inital condition
  void setVelToAnalyticalVel(AnalyticalConst3D<T, T>&);

  /// Set contact detection algorithm for particle-particle contact. Not yet implemented.
  void setContactDetection(ContactDetection<T, PARTICLETYPE>& contactDetection);
  /// Set contact detection algorithm for particle-particle contact. Not yet implemented.
  void setContactDetectionForPSys(ContactDetection<T, PARTICLETYPE>& contactDetection, int pSysNr);

  void print();

  void printDeep(std::string message="");

  /// console output number of particles at different material numbers mat
  void print(std::list<int> mat);

  /// console output of escape (E), capture (C) rate for material numbers mat
  void captureEscapeRate(std::list<int> mat);

  /// Console output of differential escape rate for material numbers mat (e.g. material of outlet).
  /// Initialisation to be done for t=0: globalPSum = 0, pSumOutlet = 0, diffEscapeRate = 0, maxDiffEscapeRate = 0
  /// Writes maximal escape rate in maxDiffEscapeRate and average rate between tStart and tEnd in avDiffEscapeRate
  /// Set genPartPerTimeStep to simulate a steady state when no new particles are generated
  void diffEscapeRate(std::list<int> mat, int& globalPSum, int& pSumOutlet, T& diffEscapeRate, T& maxDiffEscapeRate,
                      int iT, int iTConsole, T genPartPerTimeStep = 0);
  void diffEscapeRate(std::list<int> mat, int& globalPSum, int& pSumOutlet, T& diffEscapeRate, T& maxDiffEscapeRate,
                      int iT, int iTConsole, T genPartPerTimeStep,
                      T& avDiffEscapeRate, T latticeTimeStart, T latticeTimeEnd);

  /// Get Output of particleMovement
  /// Write the data of the particle movement into an txtFile
  void getOutput(std::string filename, int iT, T conversionFactorTime,
                 unsigned short particleProperties);

  /// Not relevant. But class must inherit from SuperStructure3D so we are forced to implement these functions.
  bool* operator()(int iCloc, int iX, int iY, int iZ, int iData) override
  {
    return nullptr;
  }
  int getDataSize() const override
  {
    return 0;
  }
  int getDataTypeSize() const override
  {
    return 0;
  }

  /// Particle-Fluid interaction for subgrid scale particles
  //  template<typename DESCRIPTOR>
  //  void particleOnFluid(SuperLattice3D<T, DESCRIPTOR>& sLattice, T eps, SuperGeometry3D<T>& sGeometry);
  //  template<typename DESCRIPTOR>
  //  void resetFluid(SuperLattice3D<T, DESCRIPTOR>& sLattice);

  /// returns the Stokes number
  template<typename DESCRIPTOR>
  T getStokes(UnitConverter<T, DESCRIPTOR>& conv, T pRho, T rad)
  {
    return pRho * std::pow(2.*rad, 2) * conv.getCharPhysVelocity() / (18.*conv.getCharPhysLength() *
           (conv.getPhysViscosity() * conv.getPhysDensity()));
  };

  friend class SuperParticleSysVtuWriter<T, PARTICLETYPE> ;
  friend class SuperParticleSysVtuWriterMag<T> ;

  enum particleProperties
    : unsigned short {position = 1, velocity = 2, radius = 4, mass = 8,
                    force = 16, storeForce = 32
                   };

protected:
  mutable OstreamManager clout;

  /// Init the SuperParticleSystem
  void init();
  /// Redistribute particles on compute nodes
  void updateParticleDistribution();
  /// Find the cuboid the particle is on.
  bool findCuboid(PARTICLETYPE<T>&, int overlap);
  bool findCuboid(PARTICLETYPE<T>&);
  /// Check if particle is still on cuboid
  // TODO overlap caste to int
  bool checkCuboid(PARTICLETYPE<T>& p, T overlap);
  bool checkCuboid(PARTICLETYPE<T>& p, T overlap, int iC);
  int countLocMaterial(int mat);

  /// Add a shadow particle to system
  void addShadowParticle(PARTICLETYPE<T> &p);

  /// The particleSystems. One per cuboid
  std::vector<ParticleSystem3D<T, PARTICLETYPE>*> _pSystems;
  /// The superGeometry
  SuperGeometry3D<T>& _superGeometry;
  /// Rank of neighbouring cuboids
  std::list<int> _rankNeighbours;
  /// Numbers of neighbouring cuboids
  std::vector<std::vector<int> > _cuboidNeighbours;
  // TODO: !attention! here T _overlap; class SuperStructure3D<T> has int _overlap
  // but superParticleSystem<T, PARTICLETYPE<T>> overlap is of type T
  T _overlap;
  /// temporary variables
  std::map<int, std::vector<double> > _send_buffer;
  std::multimap<int, PARTICLETYPE<T> > _relocate;
  std::multimap<int, PARTICLETYPE<T> > _relocateShadow;
};

// Magnetic particle type
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(IndicatorF3D<double>& ind, double mas,
    double rad, int no, int id,
    std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel,
    std::vector<double> torque, double magnetisation, int sActivity);
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(IndicatorF3D<double>& ind, double mas,
    double rad, int no, int id,
    std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel,
    std::vector<double> torque, double magnetisation, int sActivity);
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(IndicatorF3D<double>& ind,
    std::set<int>  material, double mas, double rad, int no, int id,
    std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel,
    std::vector<double> torque, double magnetisation, int sActivity);

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticlesdMomRandom();
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticles(std::vector<double> dMoment,
    std::vector<double> vel, std::vector<double> aVel,
    std::vector<double> torque, double magnetisation, int sActivity);
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::prepareAgglomerates();
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::initAggloParticles();
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::findAgglomerates(int iT, int itVtkOutputMagParticles);
template<>
bool SuperParticleSystem3D<double, MagneticParticle3D>::particleSActivityTest(int sActivity);
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::simulate(double dT, std::set<int> sActivityOfFreeParticle, bool scale) ;
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::simulateWithTwoWayCoupling_Mathias ( double dT,
                                    ForwardCouplingModel<double,MagneticParticle3D>& forwardCoupling,
                                    BackCouplingModel<double,MagneticParticle3D>& backCoupling,
                                    int material, int subSteps, bool resetExternalField, bool scale );
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::simulateWithTwoWayCoupling_Davide ( double dT,
                                    ForwardCouplingModel<double,MagneticParticle3D>& forwardCoupling,
                                    BackCouplingModel<double,MagneticParticle3D>& backCoupling,
                                    int material, int subSteps, bool resetExternalField, bool scale );


}  //namespace olb

#endif /* SUPERPARTICLESYSTEM_3D_H */
