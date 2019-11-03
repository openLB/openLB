/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Benjamin Förster
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


#ifndef VTI_READER_HH
#define VTI_READER_HH

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cassert>

#include "geometry/cuboid2D.h"
#include "vtiReader.h"
#include "communication/heuristicLoadBalancer.h"

namespace olb {

//template< typename T, typename BaseType> class SuperData3D;

/* ------------------ BaseVTIreader -----------------*/

template<typename T>
BaseVTIreader<T>::BaseVTIreader( const std::string& fName, int dim,
                                 std::string dName, const std::string class_name )
  : clout(class_name), _dim(dim), _size(0), _origin(_dim, 0), _extent(_dim, 0),
    _delta(0), _xmlReader(fName), _nCuboids(0)
{
  // Read WholeExtent (2 * _dim ints) from XML File and calculate _extent
  std::vector<int> wholeExtend = readExtent(&_xmlReader["ImageData"], "WholeExtent");
  _extent = getNbNodes(wholeExtend);

  // Read _delta
  std::stringstream stream_val_1(_xmlReader["ImageData"].getAttribute("Spacing"));
  stream_val_1 >> _delta;

  // Read _origin
  std::stringstream stream_val_2(_xmlReader["ImageData"].getAttribute("Origin"));
  for (auto& origin_i : _origin) {
    stream_val_2 >> origin_i;
  }

  // Read _size and count cuboids (_nCuboids)
  for ( auto& piece : _xmlReader["ImageData"] ) {
    if (piece->getName() == "Piece") {
      // Read _size from first dName PointData Tag
      if ( _nCuboids == 0 ) {
        for (auto &dataArray : (*piece)["PointData"]) {
          if (dataArray->getAttribute("Name") == dName && dataArray->getName() == "DataArray") {
            _size = this->getSize(*dataArray);
          }
        }
      }

      _nCuboids++;
    }
  }
}

template<typename T>
void BaseVTIreader<T>::printInfo()
{
  clout << "Information on VTIreader Data:" << std::endl;
  clout << "Origin: ";
  for (auto& origin_i : _origin) {
    clout << origin_i << " ";
  }
  clout << std::endl;

  clout << "Extend: ";
  for (auto& extend_i : _extent) {
    clout << extend_i << " ";
  }
  clout << std::endl;

  clout << "Spacing: " << _delta << std::endl;
}

template<typename T>
std::vector<int> BaseVTIreader<T>::readExtent(const XMLreader* reader, std::string extAttrName)
{
  // An extent is in the form of four (or six) integers "x0 x1 y0 y1 z0 z1"
  //  each representing a node number
  std::stringstream extstr(reader->getAttribute(extAttrName));
  std::vector<int> extents;
  int tmp;
  for (int i = 0; i < 2 * _dim; ++i) {
    extstr >> tmp;
    extents.push_back(tmp);
  }
  return extents;
}

template<typename T>
std::vector<int> BaseVTIreader<T>::getNbNodes(std::vector<int>& extents)
{
  // Convert 4D (or 6D) extents vector into 2D (3D) extent vector
  std::vector<int> nNodes;
  for ( int i = 0; i < _dim; i++ ) {
    nNodes.push_back(extents[ 2*i + 1 ] - extents[ 2* i ] + 1);
  }
  return nNodes;
}

template<typename T>
int BaseVTIreader<T>::getSize(const XMLreader& tag)
{
  // read the NumberOfComponents-Attribute (VTI standard) and return as integer
  int size = std::atoi((tag.getAttribute("NumberOfComponents")).c_str());
  if (size == 0) {
    clout << "VTI READER: NumberOfComponents zero or not given!" << std::endl;
    clout << "EXAMPLE: <DataArray Name='physVelocity' NumberOfComponents='3'" << std::endl;
    exit(1);
  }
  return size;
}


/* -------------------- BaseVTIreader3D --------------------- */

template<typename T, typename BaseType>
BaseVTIreader3D<T,BaseType>::BaseVTIreader3D( const std::string& fName, std::string dName,
    const std::string class_name)
  : BaseVTIreader<T>(fName, 3, dName, class_name)
{
}

template<typename T, typename BaseType>
void BaseVTIreader3D<T,BaseType>::readCuboid(Cuboid3D<T>& cuboid, XMLreader* piece)
{
  if (piece->getName() == "Piece") {
    std::vector<int> extents = this->readExtent(piece, "Extent");
    std::vector<int> extent = this->getNbNodes(extents);
    // int extents[i] is node number => multiply with _delta to get coordinate
    cuboid.init(extents[0] * this->_delta,
                extents[2] * this->_delta,
                extents[4] * this->_delta,
                this->_delta,
                extent[0],
                extent[1],
                extent[2]);
  }
}

template<typename T, typename BaseType>
bool BaseVTIreader3D<T,BaseType>::readBlockData(BlockData3D<T,BaseType>& blockData,
    const XMLreader& pieceTag, const std::string dName)
{
  /*** This is the main data reader method. All other methods use this one for accessing data ***/

  // Calculate number of cells in each direction
  std::vector<int> extents = this->readExtent(&pieceTag, "Extent");
  std::vector<int> extent = this->getNbNodes(extents);

  // Iterate through all <DataArray> tags and take the one with the given Name attribute
  for (auto & dataArray : pieceTag["PointData"]) {
    if (dataArray->getAttribute("Name") == dName && dataArray->getName() == "DataArray") {
      std::string data_str;
      if (dataArray->read(data_str)) {
        std::stringstream stream_val(data_str);

        // Careful: respect ordering in VTI File
        for (int iz = 0; iz < blockData.getNz(); iz++) {
          for (int iy = 0; iy < blockData.getNy(); iy++) {
            for (int ix = 0; ix < blockData.getNx(); ix++) {
              for (int iSize=0; iSize < blockData.getSize(); iSize++) {
                BaseType tmp;
                stream_val >> tmp;
                // write tmp into blockData
                blockData.get(ix, iy, iz, iSize) = tmp;
              }
            }
          }
        }
      }

      // set correct blockData Name (GenericF)
      //      blockData.getName() = dName;

      return true;
    }
  }
  return false;
}


/* ---------------- BlockVTIreader3D -------------------*/

template<typename T,typename BaseType>
BlockVTIreader3D<T,BaseType>::BlockVTIreader3D(const std::string& fName, const std::string dName )
  : BaseVTIreader3D<T,BaseType>(fName, dName, "BlockVTIreader3D"),
    _cuboid(this->_origin, this->_delta, this->_extent),
    _blockData(_cuboid, this->_size)
{
  // Only read the first <Piece> tag in the XML file
  this->readBlockData(_blockData, this->_xmlReader["ImageData"]["Piece"], dName);
}


template<typename T,typename BaseType>
BlockData3D<T,BaseType>& BlockVTIreader3D<T,BaseType>::getBlockData()
{
  return _blockData;
}

template<typename T,typename BaseType>
Cuboid3D<T>& BlockVTIreader3D<T,BaseType>::getCuboid()
{
  return _cuboid;
}


/* --------------- SuperVTIreader3D ---------------------*/

template<typename T,typename BaseType>
SuperVTIreader3D<T,BaseType>::~SuperVTIreader3D()
{
  delete _loadBalancer;
  delete _cGeometry;
  delete _superData;
}

template<typename T,typename BaseType>
SuperVTIreader3D<T,BaseType>::SuperVTIreader3D(const std::string& fName, const std::string dName )
  : BaseVTIreader3D<T,BaseType>(fName, dName, "SuperVTIreader3D")
{
  this->clout << "Start reading \"" << fName << "\"... "
              << "(" << this->_nCuboids << " cuboids)" << std::endl;

  // Create CuboidGeometry
  _cGeometry = new CuboidGeometry3D<T> (this->_origin, this->_delta, this->_extent);

  this->clout << "* Reading Cuboid Geometry..." << std::endl;

  // Fill CuboidGeometry
  readCuboidGeometry();

  _cGeometry->printExtended();

  // Create LoadBalancer
  _loadBalancer = new HeuristicLoadBalancer<T> (*_cGeometry);

  // Create SuperData (allocation of the data, this->_size is already known!)
  _superData = new SuperData3D<T,BaseType> ( *_cGeometry, *_loadBalancer, 2, this->_size);

  this->clout << "* Reading BlockData..." << std::endl;

  // Fill data objects
  readSuperData(dName);

  this->clout << "VTI Reader finished." << std::endl;
}

template<typename T,typename BaseType>
void SuperVTIreader3D<T,BaseType>::readCuboidGeometry()
{
  //_cGeometry->clearCuboids();
  for ( auto& piece : this->_xmlReader["ImageData"] ) {
    if (piece->getName() == "Piece") {
      std::vector<int> extents = this->readExtent(piece, "Extent");
      std::vector<int> extent = this->getNbNodes(extents);
      // int extent[i] is node number => multiply with _delta to get coordinate
      Cuboid3D<T> cuboid(extents[0] * this->_delta,
                         extents[2] * this->_delta,
                         extents[4] * this->_delta,
                         this->_delta,
                         extent[0],
                         extent[1],
                         extent[2]);
      _cGeometry->add(cuboid);
    }
  }
}

template<typename T,typename BaseType>
void SuperVTIreader3D<T,BaseType>::readSuperData(const std::string dName)
{
  int counter = 0;
  // Iterate over all <Piece> tags
  for (auto & piece : this->_xmlReader["ImageData"]) {
    if (piece->getName() == "Piece") {
      this->readBlockData(_superData->get(counter), *piece, dName);
      counter++;
    }
  }
}

template<typename T,typename BaseType>
SuperData3D<T,BaseType>& SuperVTIreader3D<T,BaseType>::getSuperData()
{
  return *_superData;
}

template<typename T,typename BaseType>
CuboidGeometry3D<T>& SuperVTIreader3D<T,BaseType>::getCuboidGeometry()
{
  return *_cGeometry;
}

template<typename T,typename BaseType>
LoadBalancer<T>& SuperVTIreader3D<T,BaseType>::getLoadBalancer()
{
  return *_loadBalancer;
}






// ************************************ old 2D ********************************************** //

/*
template<typename T>
void VTIreader2D<T>::printInfo()
{
  clout << "Origin: " << _x0 << " "  << _y0 << std::endl;
  clout << "Extend: " << _x  << " "  << _y  << std::endl;
  clout << "Spacing: " << _delta << std::endl;
}

template<typename T>
VTIreader2D<T>::VTIreader2D(const std::string& fName ) : XMLreader(fName)
{
  // get _x, _y, _z
  int x, y, z;
  std::stringstream stream_val_0((*this)["ImageData"].getAttribute("WholeExtent"));
  stream_val_0 >> x >> _x >> y >> _y >> z >> _z;
  _x = _x-x+1;
  _y = _y-y+1;
  _z = _z-z+1;

  // get _delta
  std::stringstream stream_val_1((*this)["ImageData"].getAttribute("Spacing"));
  stream_val_1 >> _delta;

  std::stringstream stream_val_2((*this)["ImageData"].getAttribute("Origin"));
  stream_val_2 >> _x0 >> _y0 >> _z0;
}

template<typename T>
VTIreader2D<T>::~VTIreader2D() {}

template<typename T>
void VTIreader2D<T>::getCuboid(Cuboid2D<T>& cuboid)
{
  cuboid.init(_x0, _y0, _delta, _x, _y);
}

template<typename T>
void VTIreader2D<T>::getCuboids(std::vector<Cuboid2D<T>* >& cuboids)
{
  std::vector<XMLreader*>::const_iterator it;
  int i = 0;
  for (it = (*this)["ImageData"].begin(); it != (*this)["ImageData"].end(); it++) {
    if ((*it)->getName() == "Piece") {
      ++i;
      std::stringstream extstr((*it)->getAttribute("Extent"));
      int extents[6];
      for (int i = 0; i < 6; ++i) {
        extstr >> extents[i];
      }
      Cuboid2D<T>* cuboid = new Cuboid2D<T>(extents[0], extents[2], _delta,
                                            extents[1]-extents[0]+1,
                                            extents[3]-extents[2]+1);
      cuboids.push_back(cuboid);
    }
  }
}
*/
/*
template<typename T>
void VTIreader2D<T>::getScalarMultiPieceData(std::vector<const ScalarFieldBase2D<T>* >& bases, const std::string dName)
{
  std::vector<XMLreader*>::const_iterator it;
  for (it = (*this)["ImageData"].begin(); it != (*this)["ImageData"].end(); it++) {
    if ((*it)->getName() == "Piece") {
      std::stringstream extstr((*it)->getAttribute("Extent"));
      int extents[6];
      for (int i = 0; i < 6; i++) {
        extstr >> extents[i];
      }
      std::vector<XMLreader*>::const_iterator it2;
      for (it2 = (*it)->operator[]("PointData").begin(); it2 != (*it)->operator[]("PointData").end(); it2++) {
        if ((*it2)->getAttribute("Name") == dName && (*it2)->getName() == "DataArray") {
          ScalarField2D<T>* tmp = new ScalarField2D<T>(extents[1]-extents[0]+1, extents[3]-extents[2]+1);
          tmp->construct();
          std::string data_str;
          if ((*it2)->read(data_str)) {
            std::stringstream stream_val(data_str);
            for (int iz = 0; iz < extents[5]-extents[4]+1; iz++) {
              for (int iy = 0; iy < extents[3]-extents[2]+1; iy++) {
                for (int ix = 0; ix < extents[1]-extents[0]+1; ix++) {
                  T tmp2;
                  stream_val >> tmp2;
                  tmp->get(ix, iy) = tmp2;
                }
              }
            }
          }
          bases.push_back(tmp);
        }
      }
    }
  }
}

template<typename T>
void VTIreader2D<T>::getVectorMultiPieceData(std::vector<const TensorFieldBase2D<T, 2>* >& bases, const std::string dName)
{
  std::vector<XMLreader*>::const_iterator it;
  for (it = (*this)["ImageData"].begin(); it != (*this)["ImageData"].end(); it++) {
    if ((*it)->getName() == "Piece") {
      std::stringstream extstr((*it)->getAttribute("Extent"));
      int extents[6];
      for (int i=0; i<6; i++) {
        extstr >> extents[i];
      }
      std::vector<XMLreader*>::const_iterator it2;
      for (it2 = (*it)->operator[]("PointData").begin(); it2 != (*it)->operator[]("PointData").end(); it2++) {
        if ((*it2)->getAttribute("Name") == dName && (*it2)->getName() == "DataArray") {
          TensorField2D<T, 2>* tmp = new TensorField2D<T, 2>(extents[1]-extents[0]+1, extents[3]-extents[2]+1);
          tmp->construct();
          std::string data_str;
          if ((*it2)->read(data_str)) {
            std::stringstream stream_val(data_str);
            for (int iz = 0; iz < extents[5]-extents[4]+1; iz++) {
              for (int iy = 0; iy < extents[3]-extents[2]+1; iy++) {
                for (int ix = 0; ix < extents[1]-extents[0]+1; ix++) {
                  T tmpval;
                  for (int i=0; i < 2; i++) {
                    stream_val >> tmpval;
                    tmp->get(ix, iy)[i] = tmpval;
                  }
                  stream_val >> tmpval;
                }
              }
            }
          }
          bases.push_back(tmp);
        }
      }
    }
  }
}

template<typename T>
bool VTIreader2D<T>::getScalarData(ScalarField2D<T>* base, const std::string dName)
{
  std::vector<XMLreader*>::const_iterator it;
  for (it = (*this)["ImageData"]["Piece"]["PointData"].begin(); it != (*this)["ImageData"]["Piece"]["PointData"].end(); it++) {
    if ((*it)->getAttribute("Name") == dName && (*it)->getName() == "DataArray") {
      ScalarField2D<T>* tmp = new ScalarField2D<T>(_x, _y);
      tmp->construct();
      std::string data_str;
      if ((*it)->read(data_str)) {
        std::stringstream stream_val(data_str);
        for (int iz = 0; iz < _z; iz++) {
          for (int iy = 0; iy < _y; iy++) {
            for (int ix = 0; ix < _x; ix++) {
              T tmp2;
              stream_val >> tmp2;
              tmp->get(ix, iy) = tmp2;
            }
          }
        }
      }
      base->swap(*tmp);
      delete tmp;
      return true;
    }
  }
  return false;
}

template<typename T>
bool VTIreader2D<T>::getVectorData(TensorField2D<T, 2>* base, const std::string dName)
{
  std::vector<XMLreader*>::const_iterator it;
  for (it = (*this)["ImageData"]["Piece"]["PointData"].begin(); it != (*this)["ImageData"]["Piece"]["PointData"].end(); it++) {
    if ((*it)->getAttribute("Name") == dName && (*it)->getName() == "DataArray") {
      TensorField2D<T, 2>* tmp = new TensorField2D<T, 2>(_x, _y);
      tmp->construct();
      std::string data_str;
      if ((*it)->read(data_str)) {
        std::stringstream stream_val(data_str);
        for (int iz = 0; iz < _z; iz++) {
          for (int iy = 0; iy < _y; iy++) {
            for (int ix = 0; ix < _x; ix++) {
              T tmpval;
              for (int i = 0; i < 2; i++) {
                stream_val >> tmpval;
                tmp->get(ix, iy)[i] = tmpval;
              }
              stream_val >> tmpval;
            }
          }
        }
      }
      base->swap(*tmp);
      delete tmp;
      return true;
    }
  }
  return false;
}
*/


}
#endif
