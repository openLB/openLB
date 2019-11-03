/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011 Lukas Baron, Mathias Krause
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


#ifndef OSTREAM_MANAGER_H
#define OSTREAM_MANAGER_H

#include <iostream>
#include <sstream>
#include <iomanip>

namespace olb {

/// userdefined stream buffer for OstreamManager
/**that prefixes each line with a user specified string in squared brackets*/
class OMBuf : public std::stringbuf {
private:
  std::ostream* output;
  std::string text;
  static bool multiOutput;
public:
  OMBuf();
  ~OMBuf() override;
  OMBuf(const OMBuf& rhs);
  OMBuf& operator=(const OMBuf& rhs);
  void swap(OMBuf& rhs);
  OMBuf(std::ostream& str, std::string classname);
  void setMultiOutput(bool b);
  /// sync the stream with the output:
  /** 1) first Output "[text] ", then the buffer,<br>
    * 2) reset the buffer<br>
    * 3) and flush the actual output stream*/
  int sync() override;
};

/// class for marking output with some text
/** The principle of this class consists of writing all output first in a userdefined Buffer of type OMBuf. On a flush it spits out at first the userdefined text in squared brackets and afterwards everything from the buffer.
An object of this class can be used (almost) exactly like a normal std::cout with the <<-Operator and std::endl.

How to implement in the code of a class with some output:
<pre>
#include<olb/io>

ExampleClass {
private:
  OstreamManager clout;
public:
  ExampleClass() : clout(std::cout,"ExampleClass")
    {}
  showOutput() {
    clout << "where am I?" << std::endl;
  }
}
</pre>
A call of the function showOutput() will show in terminal:
<pre>
[ExampleClass] Where am I?
</pre>
Please note, that the control character <b>\\n</b> - in contrast to std::endl - will not force a flush of the outstream. As consequence, the new line won't be prefixed with a user specified text. Other control character might not work either (untested).
*/
class OstreamManager : public std::ostream {
private:
  // OstreamManager clout should use it's own special buffer
  /// special, overloaded buffer
  OMBuf buffer;

public:
  // standard constructor (should be unnecessary due to absence of classname-text)
  // OstreamManager();
  /// constructor that uses std::cout by default
  OstreamManager(std::string classname);
  /// constructor for default usage
  OstreamManager(std::ostream& str, std::string classname);
  // Copy construction
  OstreamManager(const OstreamManager& rhs);
  // Copy assignment
  OstreamManager& operator=(const OstreamManager& rhs);
  // Destructor
  ~OstreamManager() override;
  /// enable message output for all MPI processes, disabled by default
  void setMultiOutput(bool b);
};

} // namespace olb

#endif
