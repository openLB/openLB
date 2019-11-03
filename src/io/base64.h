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


#ifndef BASE64_H
#define BASE64_H

#include <iosfwd>

namespace olb {

template<typename T>
class Base64Encoder {
public:
  Base64Encoder(std::ostream& ostr_, size_t fullLength_);
  void encode(const T* data, size_t length);
private:
  void fillOverflow(const unsigned char* charData, size_t charLength, size_t& pos);
  void flushOverflow();
  void writeSize();
  void encodeBlock( const unsigned char* data);
  void encodeUnfinishedBlock( const unsigned char* data, int length);
private:
  static const char enc64[65];
private:
  std::ostream& ostr;
  size_t charFullLength;
  size_t numWritten;
  int numOverflow;
  unsigned char overflow[3];
};

template<typename T>
class Base64Decoder {
public:
  Base64Decoder(std::istream& istr_, size_t fullLength_);
  void decode(T* data, size_t length);
private:
  void flushOverflow(unsigned char* charData, size_t charLength, size_t& pos);
  unsigned char getNext();
  void decodeBlock(unsigned char* data);
private:
  static const char dec64[82];
private:
  std::istream& istr;
  size_t charFullLength;
  size_t numRead;
  int posOverflow;
  unsigned char overflow[3];
};

} // namespace olb

#endif
