/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2016 Jonas Latt, Mathias J. Krause
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

#ifndef SERIALIZER_IO_HH
#define SERIALIZER_IO_HH

#include "serializerIO.h"
#include "base64.h"
#include "core/olbDebug.h"
#include <limits>
#include <istream>
#include <ostream>
#include <fstream>
#include <sstream>

namespace olb {

void serializer2ostr(Serializer& serializer, std::ostream& ostr, bool enforceUint)
{
  serializer.resetCounter();
  // write binary size into first integer of stream
  size_t binarySize = serializer.getSize();
  if (enforceUint) {
    Base64Encoder<unsigned int> sizeEncoder(ostr, 1);
    OLB_PRECONDITION(binarySize <= std::numeric_limits<unsigned int>::max());
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
  } else {
    Base64Encoder<size_t> sizeEncoder(ostr, 1);
    sizeEncoder.encode(&binarySize, 1);
  }

  Base64Encoder<bool> dataEncoder (ostr, binarySize);

  size_t blockSize;
  const bool* dataBuffer = nullptr;
  while (dataBuffer = serializer.getNextBlock(blockSize, false), dataBuffer != nullptr) {
    dataEncoder.encode(dataBuffer, blockSize);
  }
  serializer.resetCounter();
}

void istr2serializer(Serializer& serializer, std::istream& istr, bool enforceUint)
{
  //size_t binarySize = serializer.getSize();
  serializer.resetCounter();

  // read binary size from first integer of stream
  size_t binarySize;
  if (enforceUint) {
    unsigned int uintBinarySize;
    Base64Decoder<unsigned int> sizeDecoder(istr, 1);
    sizeDecoder.decode(&uintBinarySize, 1);
    binarySize = uintBinarySize;
  } else {
    Base64Decoder<size_t> sizeDecoder(istr, 1);
    sizeDecoder.decode(&binarySize, 1);
  }
  //OLB_PRECONDITION(binarySize == serializer.getSize());


  Base64Decoder<bool> dataDecoder(istr, binarySize);

  size_t blockSize;
  bool* dataBuffer = nullptr;
  while (dataBuffer = serializer.getNextBlock(blockSize, true), dataBuffer != nullptr) {
    dataDecoder.decode(dataBuffer, blockSize);
  }
  serializer.resetCounter();
}

} // namespace olb

#endif
