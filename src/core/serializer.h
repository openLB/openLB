/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Mathias J. Krause, Benjamin Förster
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


#ifndef SERIALIZER_H
#define SERIALIZER_H

#include <iostream>
#include <vector>
#include <map>
#include <utility>

namespace olb {

class Serializable;

/// Class for writing, reading, sending and receiving `Serializable` objects.
/**
 * __For detailed information on the serialization concept, see the `Serializable` documentation.__
 */
class Serializer {

private:
  /// Object to be (de-)serialized
  Serializable& _serializable;
  /// Counter for the current block number being processed
  std::size_t _iBlock;
  /// Total memory size in bits (computed by `computeSizes()`)
  std::size_t _size;
  /// Default file name for IO
  std::string _fileName;

public:
  /// Constructor
  /**
   * If `serializable` is omitted, it has to be provided in the save method.
   */
  Serializer(Serializable& serializable, std::string fileName = "");

  /// Resets the `_iBlock` counter
  void resetCounter();

  /// Returns the total memory size in bits
  std::size_t getSize() const;

  /// Returns pointer to the memory of the current block and increments `iBlock`
  bool* getNextBlock(std::size_t& sizeBlock, const bool loadingMode);

  /// Loads a file and pushes the data into the serialized class. Always in parallel, i.e. one file per rank.
  /**
   * \todo implement similar methods for sending and receiving data through MPI
   */
  bool load(std::string fileName = "", const bool enforceUint=false);

  /// Save `_serializable` into file `filename`. Always in parallel, i.e. one file per rank.
  bool save(std::string fileName = "", const bool enforceUint=false);

  /// computes `_size` based on the individual definition of `getBlock()`
  void computeSize(const bool enforceRecompute=false);

private:
  /// Set `fileName` to `_fileName` if empty and set it to `"Serializable"` if both equal ""
  void validateFileName(std::string &fileName);
  /// Returns full file name for `_fileName`
  const std::string getFullFileName(const std::string& fileName);
};


/// Base class for serializable objects of _constant size_. For _dynamic size_ use `BufferSerializable`.
/**
 * All serializable classes have to implement their individual `getBlock()` method.
 * An individual `getNblock()` method must also be provided. An individual `getSerializableSize()` method should
 * also be provided for efficiency reasons.
 *
 * The `sumNblock` and `sumSerializableSize` operator structs can be used for accumulation of `getNblock()` methods
 * (and `getSerializableSize()` respectively), e.g. with an array or a `std::vector<Serializable>`.
 *
 * All `Serializable` subclasses with _dynamic size_ (unknown at compile time, e.g. holding `std::vector` or `std::map`
 * members) have to inherit from `BufferSerializable`. _Note: If the dynamic size is computable through __constant__
 * values (see `BlockLattice2D`), the `Serializable` does not need to be a `BufferSerializable`.
 *
 * ### The Basic Serialization Concept ###
 *
 * Any serializable class inherits from either `Serializable` or `BufferSerializable` (if it contains _dynamic-sized_
 * member variables) and has to implement its individual `getBlock()` method.
 *
 * `getBlock()` is called by the `Serializer` repeatedly with an increasing counter `iBlock`. `getBlock()` returns a
 * `(bool*)` address to the _i-th_ data block and fills `sizeBlock` with the corresponding size. As long as `getBlock()`
 * does not return a `nullptr`, `iBlock` is increased and `getBlock()` is called again by `Serializer`.
 *
 * It is _strongly recommended_ (and __obligatory__ for the correct usage of `register` methods) to define
 * `std::size_t currentBlock = 0;` within `getBlock()`. `currentBlock` will be increased by the `register` methods by
 * the number of blocks they occupy.
 *
 * For user's convenience the `Serializable` class provides `register` methods for _data of constant size_:
 * Method                               | Suitable for
 * ------------------------------------ | --------------------------
 * `registerVar()`                      | Primitive data types and arrays of those (e.g. `int`, `double`, `std::string`, ...)
 * `registerSerializableOfConstSize()`  | Constant-sized `Serializable` object
 * `registerSerializablesOfConstSize()` | Array of constant-sized `Serializable` objects
 *
 * - In `registerVar()`, `currentBlock` is counted up by 1 - `arrays` are treated as one block of size
 * `sizeof(DataType) * arrayLength`.
 * - In `registerSerializableOfConstSize()`, `currentBlock` is increased by
 * `getNblock()` of the given `Serializable`.
 * - In `registerSerializablesOfConstSize()`, `currentBlock` is increased by
 * `arrayLength * data.getNblock()` of the given `Serializable`.
 *
 *
 * __Note:__ Dynamic-sized objects need to inherit from `BufferSerializable`, which uses buffers and provides
 * additional `register` methods for:
 *
 * Method                                          | Suitable for
 * ----------------------------------------------- | --------------------------
 * `registerSerializable()`                        | Dynamic-sized `Serializable` object (for constant-sized `Serializable`
 *                                                 | use `registerSerializableOfConstSize()`
 * `registerStdVectorOfVars()`                     | `std::vector<DataType>` (for primitive `DataType`, e.g. `int`,
 *                                                 | `double`, ...)
 * `registerStdVectorOfSerializablesOfConstSize()` | `std::vector<DataType>` (for `Serializable`s of constant size)
 * `registerStdVectorOfSerializables()`            | `std::vector<DataType>` (for `Serializable`s of dynamic size)
 * `registerMap()`                                 | `std::map<DataTypeKey, DataTypeValue>` (for primitive types)
 *
 */
class Serializable {
public:
  /// Returns the address of the i-th block and its size.
  /**
   * \param iBlock      Index of the block to be returned
   * \param sizeBlock   Reference to the size of the returned block
   * \return            Pointer to the current block
   *
   * Each `getBlock()` method should look like this:
   *
   *     std::size_t currentBlock = 0;
   *     bool* dataPtr = nullptr;
   *
   *     // ... register methods...
   *
   *     return dataPtr;
   */
  virtual bool* getBlock(const std::size_t iBlock, std::size_t& sizeBlock,
                         const bool loadingMode = false) = 0;

  /// Returns the number of blocks.
  /**
   * All `Serializable` classes have to implement this method.
   */
  virtual std::size_t getNblock() const = 0;

  /// Returns the binary size of the data to be saved
  /**
   * _This method must be overloaded by all child classes._
   */
  virtual std::size_t getSerializableSize() const = 0;

  /// Save `Serializable` into file `fileName`
  bool save(std::string fileName = "", const bool enforceUint=false);

  /// Load `Serializable` from file `fileName`
  bool load(std::string fileName = "", const bool enforceUint=false);

  /// Sum functor for `getNblock()` of `std::vector<Serializable>` (for `std::accumulate`)
  /**
   * _Usage_: If you have `std::vector<Serializable> v`, then use accumulate as follows:
   *
   *     std::accumulate(v.begin(), v.end(), size_t(0), Serializable::sumNblock());
   */
  struct sumNblock : public std::binary_function<size_t, Serializable, size_t> {
    std::size_t operator()(std::size_t sum, const Serializable& s)
    {
      return sum + s.getNblock();
    }
  };

  /// Sum functor for `getSerializableSize()` of `std::vector<Serializable>` (for `std::accumulate`)
  /**
   * _Usage_: If you have `std::vector<Serializable> v`, then use accumulate as follows:
   *
   *     std::accumulate(v.begin(), v.end(), size_t(0), Serializable::sumSerializableSize());
   */
  struct sumSerializableSize : public std::binary_function<size_t, Serializable, size_t> {
    std::size_t operator()(std::size_t sum, const Serializable& s)
    {
      return sum + s.getSerializableSize();
    }
  };

protected:
  /// Register _primitive data types_ (`int`, `double`, ...) or arrays of those
  /**
   * This method is suitable for all _primitive data types_ or arrays of those. The address of the data is returned
   * in combination with the size `sizeof(DataType) * arrayLength`.
   *
   * \param iBlock         `iBlock` from `getBlock()` - to determine if this is the current block.
   * \param sizeBlock      `sizeBlock` from `getBlock()` - will be filled if this is the current block.
   * \param currentBlock   _local_ variable of `getBlock()` - will always be counted up by the
   *                       number of blocks this method registers (_always 1 in this case_).
   * \param dataPtr        `dataPtr` from `getBlock()` - will be filled with pointer to the data at `iBlock` if this is
   *                       the current block.
   * \param data           Reference to the data to be registered by this method. Fills `dataPtr` with a
   *                       `(bool*)`-casted pointer to `data` if this is the current block.
   * \param arrayLength    Number of elements of `DataType` in `data`, if `data` is an array. _Defaults to 1 for
   *                       single values._
   * \tparam DataType      Type of `data`
   */
  template<typename DataType>
  void registerVar(const std::size_t iBlock, std::size_t &sizeBlock, std::size_t &currentBlock, bool *&dataPtr,
                   const DataType &data, const size_t arrayLength = 1) const
  {
    if (iBlock == currentBlock) {
      sizeBlock = sizeof(DataType) * arrayLength;
      dataPtr = (bool *) (&data);
    }
    currentBlock++;
  }


  /// Register `Serializable` object of _constant size_.
  /**
   * This method is suitable for all `Serializable` objects that are of __constant size__.
   *
   * _For information about the parameters of this method, see `registerVar()`._
   *
   * This method registers a `Serializable` by simply delegating the `getBlock()` call to the `Serializable`.
   *
   * Since those `Serializable` objects __must(!)__ have constant return value from `Serializable.getNblock()`,
   * the number of blocks is known both in reading and writing mode.
   */
  template<typename DataType>
  void registerSerializableOfConstSize(const std::size_t iBlock, std::size_t &sizeBlock, std::size_t &currentBlock,
                                       bool *&dataPtr, DataType &data, const bool loadingMode=false)
  {
    static_assert(std::is_base_of<Serializable, DataType>::value, "DataType must be a Serializable.");

    if (iBlock >= currentBlock && iBlock < currentBlock + data.getNblock()) {
      dataPtr = data.getBlock(iBlock - currentBlock, sizeBlock, loadingMode);
    }
    currentBlock += data.getNblock();
  }

  /// Register an array of `Serializable` objects of _constant size_.
  /**
   * This method is suitable for all `Serializable` objects that are of __constant size__.
   *
   * _For information about the parameters of this method, see `registerVar()`._
   *
   * This method registers an array of `Serializable`s by delegating the `getBlock()` call
   * to the corresponding `Serializable`.
   *
   * Since those `Serializable` objects __must(!)__ have constant return value from `Serializable.getNblock()`,
   * the number of blocks is known both in reading and writing mode.
   */
  template<typename DataType>
  void registerSerializablesOfConstSize(const std::size_t iBlock, std::size_t &sizeBlock, std::size_t &currentBlock,
                                        bool *&dataPtr, DataType* data, const size_t arrayLength,
                                        const bool loadingMode=false)
  {
    static_assert(std::is_base_of<Serializable, DataType>::value, "DataType must be a Serializable.");

    if ( arrayLength > 0 ) {
      if (iBlock >= currentBlock && iBlock < currentBlock + arrayLength * data[0].getNblock()) {
        size_t local_iBlock = iBlock - currentBlock;
        size_t dataBlockCount = data[0].getNblock();
        dataPtr = data[local_iBlock / dataBlockCount].getBlock(local_iBlock % dataBlockCount, sizeBlock, loadingMode);
      }
      currentBlock += arrayLength * data[0].getNblock();
    }
  }
};


/// Base class for serializable objects of _dynamic size_
/**
 * All `Serializable` subclasses with _dynamic size_ (unknown at compile time , e.g. holding `std::vector` or `std::map`
 * members) have to inherit from `BufferSerializable`
 *
 * __For detailed information on the serialization concept, see the `Serializable` documentation.__
 *
 * `BufferSerializable` provides `register` methods for:
 *
 * Method                                          | Suitable for
 * ----------------------------------------------- | --------------------------
 * `registerSerializable()`                        | Dynamic-sized `Serializable` object (for constant-sized `Serializable`
 *                                                 | use `registerSerializableOfConstSize()`
 * `registerStdVectorOfVars()`                     | `std::vector<DataType>` (for primitive `DataType`, e.g. `int`,
 *                                                 | `double`, ...)
 * `registerStdVectorOfSerializablesOfConstSize()` | `std::vector<DataType>` (for `Serializable`s of constant size)
 * `registerStdVectorOfSerializables()`            | `std::vector<DataType>` (for `Serializable`s of dynamic size)
 * `registerMap()`                                 | `std::map<DataTypeKey, DataTypeValue>` (for primitive types)
 *
 */
class BufferSerializable : public Serializable {
protected:
  /// Data buffer for data that has to be buffered between two `getBlock()` iterations.
  /**
   * This is currently only used within `registerMap()`.
   *
   * _This is only used for loading._
   */
  mutable std::vector<bool*> _dataBuffer;
  /// `std::vector` of integer buffers (e.g. for `std::vector` size) to be buffered for the whole iteration process
  /**
   * Each register method for dynamic-sized objects (e.g. `std::vector` or `std::map`) uses the size buffer to
   * provide to correctly increase the `currentBlock` variable for all following register methods.
   */
  mutable std::vector<size_t> _sizeBuffer;


  /// Register `Serializable` object of _dynamic size_.
  /**
   * This method is suitable for all `Serializable` objects of _ dynamic size_, e.g. where `getNblock()` cannot be
   * evaluated correctly on an empty object.
   *
   * _For more information about the parameters of this method, see `registerVar()`._
   *
   * This method registers a `Serializable` by buffering the dynamic `getNblock()` value and afterwards
   * delegating the `getBlock()` call to the `Serializable`.
   *
   * _Note:_ `getNblock()` only works for writing mode, since it is dynamic.
   */
  template<typename DataType>
  void registerSerializable(const std::size_t iBlock, std::size_t &sizeBlock, std::size_t &currentBlock,
                            size_t &sizeBufferIndex, bool *&dataPtr, DataType &data,
                            const bool loadingMode=false)
  {
    static_assert(std::is_base_of<Serializable, DataType>::value, "DataType must be a Serializable.");

    size_t dataBlockCount = 0;

    // hold getNblock() in sizeBuffer
    if (loadingMode) { // loading -> set to 0 and wait for reading next round
      dataBlockCount = addSizeToBuffer(iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, 0);
    } else { // saving -> save getNblock from data object
      dataBlockCount = addSizeToBuffer(iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, data.getNblock());
    }

    if (iBlock >= currentBlock) {
      if (iBlock < currentBlock + dataBlockCount) {
        dataPtr = data.getBlock(iBlock - currentBlock, sizeBlock, loadingMode);
      }
    }
  }

  /// Method for registering a `std::vector<DataType>` of primitive `DataType` (`int`, `double`, ...)
  /**
   * This method registers a vector of a primitive `DataType`. The first block holds a `size_t sizeOfVector` with the
   * size of the vector, followed by `sizeOfVector` many blocks of `DataType`.
   *
   * The total number of blocks occupied by this method is `1 + sizeOfVector`.
   *
   * \param sizeBufferIndex Index counter for size buffer. Is increased by one by this method, and the size of
   *                       the registered `std::vector` is stored in the corresponding `sizeBuffer` element.
   *
   * \param data           Reference to the data to be registered by this method.
   *                       Fills `dataPtr` with a
   *                       `(bool*)`-casted pointer (_if this is the current block_) to:
   *                       - _First block_ - number of elements in the vector
   *                       - _second block_ to _last block_ - pointer to `i-th` vector element
   *
   * _For information about the other parameters of this method, see `registerVar()` documentation._
   */
  template<typename DataType>
  void registerStdVectorOfVars(const std::size_t iBlock, std::size_t &sizeBlock, std::size_t &currentBlock,
                               size_t &sizeBufferIndex, bool *&dataPtr, std::vector<DataType> &data,
                               const bool loadingMode = false)
  {
    if (iBlock >= currentBlock) {
      // process length of data vector
      size_t sizeOfVector = addSizeToBuffer(iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, data.size());

      // resize data vector from buffer (only for loading)
      if (iBlock == currentBlock && loadingMode) {
        data.resize(_sizeBuffer[sizeBufferIndex - 1]);
      }

      if (iBlock >= currentBlock && iBlock < currentBlock + sizeOfVector) {
        sizeBlock = sizeof(DataType);
        dataPtr = (bool *) (&data[iBlock - currentBlock]);
      }
      currentBlock += sizeOfVector;
    }
  }


  /// Method for registering a `std::vector<DataType>` of constant-sized `Serializable`
  /**
   * This method registers a vector of a constant-sized `Serializable`. The first block holds a `size_t sizeOfVector` with
   * the size of the vector, the second holding the number of blocks in one `DataType`,
   * followed by `sizeOfVector` blocks of `DataType`, each of which consists of `DataType.getNblock()` subblocks.
   *
   * The `sizeBuffer` is used to store the _constant(!)_ number of subblocks (`DataType.getNblock()`) of the
   * `Serializable`.
   *
   * The total number of blocks occupied by this method is `2 + sizeOfVector*nSubBlock`.
   *
   * \param sizeBufferIndex Index counter for size buffer. Is increased by __two__ by this method, and the size of
   *                       the registered `std::vector` is stored in the corresponding _first_ `sizeBuffer` element,
   *                       the fixed number of blocks of `DataType` is stored in the _second_ one.
   *
   * \param data           Reference to the data to be registered by this method. Fills `dataPtr` with a
   *                       `(bool*)`-casted pointer (_if this is the current block_) to:
   *                       - _First block_ - number of elements in the vector
   *                       - _second block_ - fixed number of block in `DataType`
   *                       - _third block_ to _last block_ - pointer to `i-th` vector element
   *
   *
   * _For information about the other parameters of this method, see `registerVar()` documentation._
   */
  template<typename DataType>
  void registerStdVectorOfSerializablesOfConstSize(const std::size_t iBlock, std::size_t &sizeBlock,
      std::size_t &currentBlock, size_t &sizeBufferIndex, bool *&dataPtr,
      std::vector<DataType> &data, const bool loadingMode = false)
  {
    static_assert(std::is_base_of<Serializable, DataType>::value, "DataType must be a Serializable.");

    if (iBlock >= currentBlock) {
      // process length of data vector
      size_t sizeOfVector = addSizeToBuffer(iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, data.size());

      // resize data vector from buffer (only for loading)
      if (iBlock == currentBlock && loadingMode) {
        data.resize(_sizeBuffer[sizeBufferIndex - 1]);
      }


      // process Serializables
      if (iBlock >= currentBlock && sizeOfVector > 0) {
        for ( DataType& dataValue : data ) {
          registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, dataValue, loadingMode);
        }
      }
    }
  }



  /// Method for registering a `std::vector<DataType>` of dynamic-sized `DataType`
  /**
   * This method registers a vector of a dynamic-sized `Serializable`. The first block holds a `size_t sizeOfVector` with
   * the size of the vector, followed by `sizeOfVector` many `registerSerializable()` calls, which occupy
   * `data[i].getNblock()` blocks each. `getNblock()` may differ for any `data[i]`.
   *
   * The `sizeBuffer` is used to store the length of the registered `std::vector` as well as the number of blocks
   * occupied by each `Serializable` in the vector. This method occupies `1 + sizeOfVector` many `sizeBuffer` elements.
   *
   * The total number of blocks occupied by this method is `1 + sizeOfVector * (1 + data[i].getNblock())`.
   *
   * \param sizeBufferIndex Index counter for size buffer. Is increased by __two__ by this method, and the size of
   *                       the registered `std::vector` is stored in the corresponding _first_ `sizeBuffer` element,
   *                       the fixed number of blocks of `DataType` is stored in the _second_ one.
   *
   * \param data           Reference to the data to be registered by this method. Fills `dataPtr` with a
   *                       `(bool*)`-casted pointer (_if this is the current block_) to:
   *                       - _First block_ - number of elements in the vector
   *                       - For each `std::vector` element:
   *                           - _First Subblock_ - number of blocks occupied by this `Serializable` object
   *                           - _second block_ to _last block_ - pointer to `i-th` block of `data[i]`
   *
   *
   * _For information about the other parameters of this method, see `registerVar()` documentation._
   */
  template<typename DataType>
  void registerStdVectorOfSerializables(const std::size_t iBlock, std::size_t &sizeBlock, std::size_t &currentBlock,
                                        size_t &sizeBufferIndex, bool *&dataPtr, std::vector<DataType> &data,
                                        const bool loadingMode = false)
  {
    static_assert(std::is_base_of<Serializable, DataType>::value, "DataType must be a Serializable.");

    if (iBlock >= currentBlock) {
      // process length of data vector
      size_t sizeOfVector = addSizeToBuffer(iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, data.size());

      // resize data vector from buffer (only for loading)
      if (iBlock == currentBlock && loadingMode) {
        data.resize(_sizeBuffer[sizeBufferIndex - 1]);
      }

      // process Serializables
      if (iBlock >= currentBlock && sizeOfVector > 0) {
        for ( DataType& dataValue : data ) {
          registerSerializable(iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, dataValue, loadingMode);
        }
      }
    }
  }



  /// Method for registering a `std::map<DataTypeKey, DataTypeValue>` of fixed-sized types (i.e. `int`, `double`)
  /**
   * This method registers a map of a fixed-sized `std::map` (consisting of `std::pairs<DataTypeKey, DataTypeValue>`).
   * The first block holds a `size_t sizeOfMap` with the size of the map,
   * followed by `sizeOfMap` blocks of `std::pairs`.
   *
   * In case of _loading_ the data, the `dataBuffer` has to be used in order to provide a valid pointer to a
   * newly created `std::pair` before inserting that same pair into the map in the next iteration.
   *
   * The total number of blocks occupied by this method is `1 + sizeOfMap`.
   *
   * \param sizeBufferIndex see `registerStdVector()`
   *
   * \param data           Reference to the data to be registered by this method. Fills `dataPtr` with a
   *                       `(bool*)`-casted pointer (_if this is the current block_) to:
   *                       - _First block_ - number of elements in the map
   *                       - _second block_ to _last block_ - pointer to `i-th` map element
   *
   * __Note:__ In _writing mode_, `dataPtr` holds a pointer to the `i-th` map element (which is a `std::pair`).
   * In _reading mode_, `dataPtr` holds a pointer to the bool* buffer, which holds a newly created `std::pair` to be
   * filled and that pair is inserted into `data` in the following round.
   *
   * For information about the other parameters of this method, see `registerVar()` documentation.
   */
  template<typename DataTypeKey, typename DataTypeValue>
  void registerMap(const std::size_t iBlock, std::size_t &sizeBlock, std::size_t &currentBlock, size_t &sizeBufferIndex,
                   bool *&dataPtr, std::map<DataTypeKey, DataTypeValue> &data, const bool loadingMode = false)
  {
    if (iBlock >= currentBlock) {
      // process length of data map
      size_t sizeOfMap = addSizeToBuffer(iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, data.size());

      if (iBlock >= currentBlock && iBlock < currentBlock + sizeOfMap + 1) {
        // determine size of pair
        sizeBlock = sizeof(std::pair<DataTypeKey, DataTypeValue>);

        // LOADING MODE
        if (loadingMode) {
          // If pair in dataBuffer => insert into data and delete
          if (iBlock > currentBlock) {
            std::pair<DataTypeKey, DataTypeValue> *pairPtr = (std::pair<DataTypeKey, DataTypeValue> *) _dataBuffer.back();
            data.insert(*pairPtr);  // copy pair into map
            delete pairPtr;         // delete pair object that was created (with new!) in the buffer
            _dataBuffer.pop_back(); // remove pointer to deleted pair from buffer vector
          }

          // push new pair into buffer and return pointer
          if (iBlock < currentBlock + sizeOfMap) {
            _dataBuffer.push_back((bool *) new std::pair<DataTypeKey, DataTypeValue>);
            dataPtr = _dataBuffer.back();
          }
        }

        // SAVING MODE
        else {
          if (iBlock < currentBlock + sizeOfMap) {
            // advance through iterator to n-th element and return pointer to pair
            auto map_it = data.begin();
            std::advance(map_it, iBlock - currentBlock);
            dataPtr = (bool *) (&(*map_it));
          }
        }
      }
      currentBlock += sizeOfMap;
    }
  }


  /// Add a `size_t` to the `sizeBuffer` in the `n-th` round and return that `size_t` in all successive rounds
  /**
   *   - increase `currentBlock` by one
   *   - increase `sizeBufferIndex` by one.
   *   - `n-th` round: push given size_t to sizeBuffer and provide pointer to it.
   */
  size_t addSizeToBuffer(const std::size_t iBlock, std::size_t& sizeBlock, std::size_t&currentBlock,
                         size_t& sizeBufferIndex, bool*& dataPtr, const size_t data) const
  {
    size_t returnSize = 0;

    if (iBlock == currentBlock) {
      // write size into _sizeBuffer vector
      _sizeBuffer.push_back(*new size_t(data));
    }

    if (iBlock >= currentBlock) {
      returnSize = _sizeBuffer[sizeBufferIndex];
    }

    // register size as var
    registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _sizeBuffer[sizeBufferIndex]);
    sizeBufferIndex++;

    return returnSize;
  }
};

} // namespace olb

#endif
