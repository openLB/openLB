/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010 Jonas Latt, Jonas Fietz, Mathias Krause
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
 * Input/Output in XML format -- header file.
 */
#ifndef XML_IO_H
#define XML_IO_H
#ifdef ADT
template <class T, unsigned DIM> class ADf;
#endif
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <typeinfo>
#include "external/tinyxml/tinyxml.h"
#include "io/ostreamManager.h"

namespace olb {

class XMLreader {
public:
  /**
   * Constructs a new XMLreader from another XMLreader
   * \param pParent The new root node for the XMLreader
   */
  XMLreader( TiXmlNode* pParent );
  /// Constructs a new XMLreader from a XML file fName
  XMLreader( const std::string& fName );
  /// destructor
  ~XMLreader();
  /// Prints out the XML structure read in, mostly for debugging purposes
  void print(int indent) const;
  /**
   * Read a value from the xml file
   * \param reference to return the value
   * \return returns the value
   */
  //bool read(bool& value, bool verbose = true) const;
  template <typename T> bool read(T& value, bool verboseOn = true, bool exitIfMissing=false) const;
#ifdef ADT
  template <typename T,unsigned DIM> bool read(ADf<T,DIM>& value, bool verboseOn = true, bool exitIfMissing=false) const;
#endif
  template <typename T> bool read(std::vector<T>& value, bool verboseOn = true, bool exitIfMissing=false) const;
  template <typename T> T get(bool verboseOn = true, bool exitIfMissing=false) const;
  /// \return a Subtree placed at name \param name The name from which to take the subtree
  XMLreader const& operator[] (std::string name) const;
  /**
   * Returns an iterator.begin() of the child XMLreader
   * This means an iterator to the next level on an XML tree.
   */
  std::vector<XMLreader*>::const_iterator begin() const;
  /**
   * Returns an iterator.end() of the child XMLreader
   * This means an iterator to the next level on an XML tree.
   */
  std::vector<XMLreader*>::const_iterator end() const;
  /// switch warnings on/off
  void setWarningsOn(bool warnings) const;
  /// return the name of the element
  std::string getName() const;
  /// \return the value of attribute
  std::string getAttribute(const std::string& aName) const;
  /// print warning if verbose mode is on
  void printWarning(std::string typeName, std::string value, bool verboseOn, bool exitIfMissing) const;
private:
  void mainProcessorIni(TiXmlNode* pParent);
  void slaveProcessorIni();
  XMLreader();
private:
  mutable bool _warningsOn;
  std::string _text;
  std::string _name;
  static XMLreader _notFound;
protected:
  mutable OstreamManager clout;
  std::map<std::string, std::string> _attributes;
  std::vector<XMLreader*> _children;
};

// methods with template

#ifdef ADT
template <typename T, unsigned DIM>
bool XMLreader::read(ADf<T,DIM>& value, bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  T tmp = T();
  if (!(valueStr >> tmp)) {
//    if ( _verboseOn ) {
//      clout << std::string("Error: cannot read value from XML element ") << _name << std::endl;
//    }
    printWarning("ADf vector", "", verboseOn, exitIfMissing);
    return false;
  }
  value = ADf<T,DIM>(tmp);
  return true;
}
#endif

template <typename T>
bool XMLreader::read(std::vector<T>& values, bool verboseOn, bool exitIfMissing ) const
{
  std::stringstream multiValueStr(_text);
  std::string word;
  std::vector<T> tmp(values);
  while (multiValueStr>>word) {
    std::stringstream valueStr(word);
    T value;
    if (!(valueStr >> value)) {
//      if ( verboseOn ) {
//        clout << std::string("Error: cannot read value array from XML element ") << _name << std::endl;
//      }
      printWarning("std::vector", "", verboseOn, exitIfMissing);
      return false;
    }
    tmp.push_back(value);
  }
  values.swap(tmp);
  return true;
}

template <typename T>
T XMLreader::get(bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  T tmp = T();
  if (!(valueStr >> tmp)) {
//    if ( verboseOn ) {
//      clout << "Error: cannot read value from XML element " << _name << std::endl;
//    }
    printWarning(typeid(T).name(), "", verboseOn, exitIfMissing);
  }
  return tmp;
}

}  // namespace olb

#endif  // XML_IO_H
