/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias Krause, Peter Weisbrod
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


#ifndef LOAD_BALANCER_H
#define LOAD_BALANCER_H

#include <vector>
#include <map>
#include "core/singleton.h"
#include "core/serializer.h"
#include "io/xmlReader.h"

namespace olb {


//template<typename T, typename BaseType> class Vector;
template<typename T> class CuboidGeometry3D;
template<typename T> class CuboidGeometry2D;
template<typename T> class HeuristicLoadBalancer;

/** Base class for all LoadBalancer.
 *  Sketch: assume we have 6 cuboids and 2 threads. Thread number 1 owns cuboid 0 and 1.
 *  Thread number 2 owns cuboid 2, 3, 4 and 5.
 *  Then we get the following configuration:
 *
 *  global cuboid number:               0   1   2   3   4   5
 *  local cuboid number of thread 0:    0   1
 *  local cuboid number of thread 1:            0   1   2   3
 *
 *  \param _glob  is a vector from 0,1,...,numberOfCuboids-1
 *  \param _loc   indicates local cuboid number in actual thread, for given global cuboid number
 *  \param _rank  indicates the processing thread of a global cuboid number
 *
 */
template<typename T>
class LoadBalancer : public BufferSerializable {
protected:
  /// number of cuboids after shrink -1 in appropriate thread
  int _size;
  /// maps global cuboid to (local) thread cuboid
  std::map<int,int> _loc;
  /// content is 0,1,2,...,_size
  std::vector<int>  _glob;
  /// maps global cuboid number to the processing thread
  std::map<int,int> _rank;
public:
  /// Default empty constructor
  LoadBalancer(int size=1);
  /// Load constructor
  LoadBalancer(int size, std::map<int,int>& loc, std::vector<int>& glob, std::map<int,int>& rank);
  /// Default empty destructor
  virtual ~LoadBalancer();
  /// Swap method
  void swap(LoadBalancer<T>& loadBalancer);
  /// returns whether `glob` is on this process
  bool isLocal(const int& glob);
  /// \return local cuboid number of relevant thread
  int loc(const int& glob);
  /// \return local cuboid number of relevant thread
  int loc(int glob) const;
  /// \return global cuboid number of given local cuboid
  int glob(int loc) const;
  /// \param glob is the global cuboid number \return rank that owns the given global cuboid number
  int rank(const int& glob);
  /// \param glob is the global cuboid number \return rank that owns the given global cuboid number
  int rank(int glob) const;
  /// \return read only acess to _size
  int size() const;

  /// equal operator
  bool operator==(const LoadBalancer<T>& rhs) const;

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

  void print(bool multiOutput = false) const;

  virtual void reInit(CuboidGeometry3D<T>& cGeometry3d, const double ratioFullEmpty=1., const double weightEmpty=.0) {}

  virtual void reInit(CuboidGeometry2D<T>& cGeometry2d, const double ratioFullEmpty=1., const double weightEmpty=.0) {}

  /// Write itself into Stringstream
  virtual void writeToStream(std::ostream& stream)
  {
    std::cout << "Error: LoadBalancer.writeToStream not implemented yet." << std::endl;
    exit(-1);
  };
  /// Write itself into File
  virtual void writeToFile(std::string fileName)
  {
    std::cout << "Error: LoadBalancer.writeToFile not implemented yet." << std::endl;
    exit(-1);
  };

};


/// Creator Function for LoadBalancer from XMLreader.
/// * LoadBalancer Data may be either in an extra file (given by the "file" attribute of the xmlReader)
//    or within the reader (XML tag) itself. Either choice is saved in lbXml.
/// * LoadBalancer Mode is determined by the "mode" attribute of lbXml.
template<typename T>
LoadBalancer<T>* createLoadBalancer(XMLreader const& xmlReader, CuboidGeometry3D<T>* cGeo = NULL)
{
  OstreamManager clout(std::cout, "createLoadBalancer");
  std::string defaultMode = "Block";

  LoadBalancer<T>* lb; // The LoadBalancer to be returned

  // Read file attribute to decide if a new XML reader has to be created
  XMLreader const* lbXml;
  std::string fileAttr = xmlReader.getAttribute("file");
  bool newFile = ( fileAttr != "Attribute not found." );
  if ( newFile ) {
    lbXml = new XMLreader(fileAttr);
  } else {
    lbXml = &xmlReader;
  }

  bool verbose = false;
  (*lbXml).setWarningsOn(false);

  std::string mode = lbXml->getAttribute("mode");
  if ( mode == "Attribute not found.") {
    clout << "Warning: Cannot read parameter from Xml-file: Mode. Set default: mode = " << defaultMode << std::endl;
    mode = defaultMode;
  }

  // Heuristic Mode
  if ( mode == "Heuristic" ) {
    // only read ratioFullEmpty - Heuristic LB will constructed from cuboidGeometry deterministicly
    double ratioFullEmpty;
    if (!(*lbXml)["RatioFullEmpty"].read<double>(ratioFullEmpty, verbose)) {
      lb = new HeuristicLoadBalancer<T>(*cGeo);
    } else {
      lb = new HeuristicLoadBalancer<T>(*cGeo, ratioFullEmpty);
    }
  }

  // Base Mode
  else if ( mode == "Base" ) {
    //LoadBalancer<T>* loadBalancer = new LoadBalancer<T>();
    //loadBalancer->load(fileName);
    clout << "LOADING BASE LB NOT IMPLEMENTED YET!" << std::endl;
    (*lbXml).setWarningsOn(true);
    //return loadBalancer;
    lb = NULL;
  }

  // Block Mode
  else if ( mode == "Block" ) {
    int size = 1;
    if (!(*lbXml)["Size"].read<int>(size, verbose)) {
      clout << "Warning: Cannot read parameter from Xml-file: Size. Set default: size = 1"
            << std::endl;
    }
    lb = new LoadBalancer<T>(size);
  }


  (*lbXml).setWarningsOn(true);

  // Delete XMLreader if it was created for the LB file
  if ( newFile ) {
    delete lbXml;
  }

  return lb;
}


/// Creator Function for LoadBalancer from fileName
template<typename T>
LoadBalancer<T>* createLoadBalancer(std::string const& fileName, CuboidGeometry3D<T>* cGeo = NULL)
{
  std::string fname = singleton::directories().getLogOutDir() + fileName + ".xml";
  XMLreader lbReader(fname);
  return createLoadBalancer(lbReader["LoadBalancer"], cGeo);
}

}  // namespace olb

#endif
