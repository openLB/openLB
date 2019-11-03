/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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

/* Acknowledgment: The strategy adopted here to encode
 * and decode Base64, and in particular the expression of the
 * arrays Base64Encoder::enc64 and Base64Decoder::dec64,
 * is inspired by the open source library b64 by Bob Trower,
 * which is distributed with a MIT license at the address
 * http://base64.sourceforge.net/b64.c
 */

#ifndef BASE64_HH
#define BASE64_HH

#include "base64.h"
#include "core/olbDebug.h"
#include <ostream>
#include <istream>


namespace olb {

////////////// class Base64Encoder ////////////////////////////////////////////

template<typename T>
const char Base64Encoder<T>::enc64[65]=
  "ABCDEFGHIJKLMNOPQRSTUVWXYZ" \
  "abcdefghijklmnopqrstuvwxyz0123456789+/";

template<typename T>
Base64Encoder<T>::Base64Encoder(std::ostream& ostr_, size_t fullLength_)
  : ostr(ostr_),
    charFullLength(fullLength_ * sizeof(T)),
    numWritten(0),
    numOverflow(0)
{ }

template<typename T>
void Base64Encoder<T>::encode(const T* data, size_t length)
{
  const unsigned char* charData = reinterpret_cast<const unsigned char*>(data);
  size_t charLength = length * sizeof(T);
  OLB_PRECONDITION( numWritten+charLength <= charFullLength );

  size_t pos=0;
  fillOverflow(charData, charLength, pos);
  while (pos+3 <= charLength) {
    encodeBlock(charData+pos);
    pos += 3;
  }
  fillOverflow(charData, charLength, pos);
  numWritten += charLength;
  if (numWritten == charFullLength) {
    flushOverflow();
  }
}

template<typename T>
void Base64Encoder<T>::fillOverflow(const unsigned char* charData, size_t charLength, size_t& pos)
{
  while (numOverflow < 3 && pos < charLength) {
    overflow[numOverflow] = charData[pos];
    ++numOverflow;
    ++pos;
  }
  if (numOverflow == 3) {
    encodeBlock(overflow);
    numOverflow = 0;
  }
}

template<typename T>
void Base64Encoder<T>::flushOverflow()
{
  if (numOverflow > 0) {
    for (int iOverflow = numOverflow; iOverflow<3; ++iOverflow) {
      overflow[iOverflow] = 0;
    }
    encodeUnfinishedBlock(overflow, numOverflow);
    numOverflow = 0;
  }
}

template<typename T>
void Base64Encoder<T>::encodeBlock( const unsigned char* data)
{
  ostr << enc64[ data[0] >> 2 ];
  ostr << enc64[ ((data[0] & 0x03) << 4) | ((data[1] & 0xf0) >> 4) ];
  ostr << (unsigned char) (enc64[ ((data[1] & 0x0f) << 2) | ((data[2] & 0xc0) >> 6) ]);
  ostr << (unsigned char) (enc64[ data[2] & 0x3f ]);
}

template<typename T>
void Base64Encoder<T>::encodeUnfinishedBlock( const unsigned char* data, int length )
{
  ostr << enc64[ data[0] >> 2 ];
  ostr << enc64[ ((data[0] & 0x03) << 4) | ((data[1] & 0xf0) >> 4) ];
  ostr << (unsigned char) (
         length == 2 ?
         enc64[ ((data[1] & 0x0f) << 2) | ((data[2] & 0xc0) >> 6) ] :
         '='
       );
  ostr << (unsigned char) ( '=' );
}


////////////// struct Base64Decoder ////////////////////////////////////////////

template<typename T>
const char Base64Decoder<T>::dec64[82]=
  "|###}rstuvwxyz{#######>?@"\
  "ABCDEFGHIJKLMNOPQRSTUVW######XYZ"\
  "[\\]^_`abcdefghijklmnopq";


template<typename T>
Base64Decoder<T>::Base64Decoder(std::istream& istr_, size_t fullLength_)
  : istr(istr_),
    charFullLength(fullLength_ * sizeof(T)),
    numRead(0),
    posOverflow(3)
{ }

template<typename T>
void Base64Decoder<T>::decode(T* data, size_t length)
{
  unsigned char* charData = reinterpret_cast<unsigned char*>(data);
  size_t charLength = length * sizeof(T);
  OLB_PRECONDITION( numRead+charLength <= charFullLength );

  size_t pos = 0;
  flushOverflow(charData, charLength, pos);
  while (pos+3 <= charLength) {
    decodeBlock(charData+pos);
    pos += 3;
  }
  if (pos < charLength) {
    decodeBlock(overflow);
    posOverflow=0;
    flushOverflow(charData, charLength, pos);
  }
  numRead += charLength;
}

template<typename T>
void Base64Decoder<T>::flushOverflow(unsigned char* charData, size_t charLength, size_t& pos)
{
  while (posOverflow < 3 && pos < charLength) {
    charData[pos] = overflow[posOverflow];
    ++pos;
    ++posOverflow;
  }
}

template<typename T>
unsigned char Base64Decoder<T>::getNext()
{
  unsigned char nextChar;
  istr >> nextChar;
  return (unsigned char) (dec64[nextChar - 43] - 62);
}

template<typename T>
void Base64Decoder<T>::decodeBlock(unsigned char* data)
{
  unsigned char input[4];
  input[0] = getNext();
  input[1] = getNext();
  input[2] = getNext();
  input[3] = getNext();
  data[0] = (unsigned char) (input[0] << 2 | input[1] >> 4);
  data[1] = (unsigned char) (input[1] << 4 | input[2] >> 2);
  data[2] = (unsigned char) (((input[2] << 6) & 0xc0) | input[3]);
}

}  // namespace olb

#endif


