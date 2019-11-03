/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 The OpenLB project
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
 * Wrapper functions that simplify the use of MPI
 */

#ifndef MPI_MANAGER_H
#define MPI_MANAGER_H

#ifdef PARALLEL_MODE_MPI
#include "mpi.h"
#include <vector>
#endif
#include <string>
#include "io/ostreamManager.h"
#ifdef ADT
#include "opti/aDiff.h"
#endif



namespace olb {

template <typename T, typename BaseType> class BlockData2D;

namespace singleton {


#ifdef PARALLEL_MODE_MPI

/// Helper class for non blocking MPI communication

class MpiNonBlockingHelper {
private:
  /// Size of the vector _mpiRequest/_mpiStatus
  unsigned _size;
  /// vector of MPI_Request
  MPI_Request *_mpiRequest;
  /// vector of MPI_Status
  MPI_Status *_mpiStatus;
public:
  /// Constructor
  MpiNonBlockingHelper();
  /// Copy construction
  MpiNonBlockingHelper(MpiNonBlockingHelper const& rhs);
  /// Copy assignment
  MpiNonBlockingHelper operator=(MpiNonBlockingHelper rhs);
  /// Destructor
  ~MpiNonBlockingHelper();

  /// Allocates memory
  void allocate(unsigned i);
  /// Frees memory
  void free();

  /// Returns the size of the vector _mpiRequest/_mpiStatus
  unsigned const& get_size() const;
  /// Read and write access _mpiRequest
  MPI_Request* get_mpiRequest() const;
  /// Get the specified request (in DEBUG mode includes a size check)
  MPI_Request* get_mpiRequest(int i) const;
  /// Read and write access _mpiStatus
  MPI_Status* get_mpiStatus() const;

  /// Swap method
  void swap(MpiNonBlockingHelper& rhs);
};

/// Wrapper functions that simplify the use of MPI

class MpiManager {
public:
  /// Initializes the mpi manager
  void init(int *argc, char ***argv);
  /// Returns the number of processes
  int getSize() const;
  /// Returns the process ID
  int getRank() const;
  /// Returns process ID of main processor
  int bossId() const;
  /// Tells whether current processor is main processor
  bool isMainProcessor() const;
  /// Returns universal MPI-time in seconds
  double getTime() const;

  /// Synchronizes the processes
  void barrier(MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data at *buf, blocking
  template <typename T>
  void send(T *buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void send(ADf<T,DIM> *buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#endif

  /// Sends data at *buf, non blocking
  template <typename T>
  void iSend(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void iSend(ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#endif

  /// Sends data at *buf, non blocking and buffered
  template <typename T>
  void ibSend(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void ibSend(ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#endif

  /// Sends data at *buf, non blocking and request free
  template <typename T>
  void iSendRequestFree(T *buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Receives data at *buf, blocking
  template <typename T>
  void receive(T *buf, int count, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void receive(ADf<T,DIM> *buf, int count, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#endif

  /// Receives data at *buf, non blocking
  template <typename T>
  void iRecv(T *buf, int count, int source, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void iRecv(ADf<T,DIM> *buf, int count, int source, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#endif

  /// Send and receive data between two partners
  template <typename T>
  void sendRecv(T *sendBuf, T *recvBuf, int count, int dest, int source, int tag = 0,
                MPI_Comm comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void sendRecv(ADf<T,DIM> *sendBuf, ADf<T,DIM> *recvBuf, int count, int dest, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
#endif

  /// Sends data to master processor
  template <typename T>
  void sendToMaster(T* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);

  /// Scatter data from one processor over multiple processors
  template <typename T>
  void scatterV(T *sendBuf, T *recvBuf, int* sendCounts, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Gather data from multiple processors to one processor
  template <typename T>
  void gatherV(T* sendBuf, T* recvBuf, int *recvCounts, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);


  /// Broadcast data from one processor to multiple processors
  template <typename T>
  void bCast(T* sendBuf, int sendCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void bCast(ADf<T,DIM>* sendBuf, int sendCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
#endif

  /// Broadcast data when root is unknown to other processors
  template <typename T>
  void bCastThroughMaster(T* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void bCastThroughMaster(ADf<T,DIM>* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);
#endif

  /// Special case for broadcasting strings. Memory handling is automatic.
  void bCast(std::string& message, int root = 0);
  /// Special case for broadcasting BlockData2D
  void bCast(BlockData2D<double,double>& sendData, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Reduction operation toward one processor
  template <typename T>
  void reduce(T& sendVal, T& recvVal, MPI_Op op, int root = 0, MPI_Comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void reduce(ADf<T,DIM>& sendVal, ADf<T,DIM>& recvVal, MPI_Op op, int root = 0, MPI_Comm = MPI_COMM_WORLD);
#endif

  /// Element-per-element reduction of a vector of data
  template <typename T>
  void reduceVect(std::vector<T>& sendVal, std::vector<T>& recvVal,
                  MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Reduction operation, followed by a broadcast
  template <typename T>
  void reduceAndBcast(T& reductVal, MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
#ifdef ADT
  template <typename T,unsigned DIM> void reduceAndBcast(ADf<T,DIM>& reductVal, MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
#endif

  /// Complete a non-blocking MPI operation
  void wait(MPI_Request* request, MPI_Status* status);

  /// Complete a series of non-blocking MPI operations
  void waitAll(MpiNonBlockingHelper& mpiNbHelper);

private:
  /// Implementation code for Scatter
  template <typename T>
  void scatterv_impl(T *sendBuf, int* sendCounts, int* displs,
                     T* recvBuf, int recvCount, int root, MPI_Comm comm);

  /// Implementation code for Gather
  template <typename T>
  void gatherv_impl(T* sendBuf, int sendCount, T* recvBuf, int* recvCounts, int* displs,
                    int root, MPI_Comm comm);
private:
  MpiManager();
  ~MpiManager();
private:
  int numTasks, taskId;
  bool ok;
  mutable OstreamManager clout;

  friend MpiManager& mpi();
};

#else

class MpiManager {
public:
  /// Initializes the mpi manager
  void init(int *argc, char ***argv, bool verbose=false) { }
  /// Returns the number of processes
  int getSize() const
  {
    return 1;
  }
  /// Returns the process ID
  int getRank() const
  {
    return 0;
  }
  /// Returns process ID of main processor
  int bossId() const
  {
    return 0;
  }
  /// Tells whether current processor is main processor
  bool isMainProcessor() const
  {
    return true;
  }

  /// Synchronizes the processes
  void barrier() const {};

  friend MpiManager& mpi();
};

#endif  // PARALLEL_MODE_MPI

inline MpiManager& mpi()
{
  static MpiManager instance;
  return instance;
}

}  // namespace singleton


}  // namespace olb


#endif  // MPI_MANAGER_H
