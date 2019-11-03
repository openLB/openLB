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
 * Wrapper functions that simplify the use of MPI, template instatiations
 */

#ifdef PARALLEL_MODE_MPI

#include "communication/mpiManager.h"
#include "core/blockData2D.h"
#include "core/olbDebug.h"
#include <iostream>
#include <algorithm>

#ifndef OLB_PRECOMPILED
#include "core/blockData2D.hh"
#endif


namespace olb {

namespace singleton {

MpiManager::MpiManager() : ok(false), clout(std::cout,"MpiManager")
{ }

MpiManager::~MpiManager()
{
  if (ok) {
    MPI_Finalize();
    ok = false;
  }
}

void MpiManager::init(int *argc, char ***argv)
{

  int ok1 = MPI_Init(argc, argv);
  int ok2 = MPI_Comm_rank(MPI_COMM_WORLD,&taskId);
  int ok3 = MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
  ok = (ok1==0 && ok2==0 && ok3==0);
  clout << "Sucessfully initialized, numThreads=" << getSize() << std::endl;
}

int MpiManager::getSize() const
{
  return numTasks;
}

int MpiManager::getRank() const
{
  return taskId;
}

int MpiManager::bossId() const
{
  return 0;
}

bool MpiManager::isMainProcessor() const
{
  return bossId() == getRank();
}

double MpiManager::getTime() const
{
  if (!ok) {
    return 0.;
  }
  return MPI_Wtime();
}

void MpiManager::barrier(MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Barrier(comm);
}

template <>
void MpiManager::send<bool>(bool *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm);
}

template <>
void MpiManager::send<char>(char *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm);
}

template <>
void MpiManager::send<int>(int *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm);
}

template <>
void MpiManager::send<float>(float *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm);
}

template <>
void MpiManager::send<double>(double *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm);
}

template <>
void MpiManager::iSend<bool>
(bool *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<char>
(char *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<int>
(int *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<float>
(float *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<double>
(double *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm, request);
  }
}


template <>
void MpiManager::ibSend<bool>
(bool *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::ibSend<char>
(char *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm, request);
  }
}

template <>
void MpiManager::ibSend<int>
(int *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::ibSend<float>
(float *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::ibSend<double>
(double *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm, request);
  }
}


template <>
void MpiManager::iSendRequestFree<bool>
(bool *buf, int count, int dest, int tag, MPI_Comm comm)
{
  MPI_Request request;
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm, &request);
  }
  MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<char>
(char *buf, int count, int dest, int tag, MPI_Comm comm)
{
  MPI_Request request;
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm, &request);
  }
  MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<int>
(int *buf, int count, int dest, int tag, MPI_Comm comm)
{
  MPI_Request request;
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm, &request);
  }
  MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<float>
(float *buf, int count, int dest, int tag, MPI_Comm comm)
{
  MPI_Request request;
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm, &request);
  }
  MPI_Request_free(&request);
}

template <>
void MpiManager::iSendRequestFree<double>
(double *buf, int count, int dest, int tag, MPI_Comm comm)
{
  MPI_Request request;
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm, &request);
  }
  MPI_Request_free(&request);
}

template <>
void MpiManager::receive<bool>(bool *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_BYTE, source, tag, comm, &status);
}


template <>
void MpiManager::receive<char>(char *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, comm, &status);
}

template <>
void MpiManager::receive<int>(int *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_INT, source, tag, comm, &status);
}

template <>
void MpiManager::receive<float>(float *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, comm, &status);
}

template <>
void MpiManager::receive<double>(double *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, comm, &status);
}

template <>
void MpiManager::sendToMaster<bool>(bool* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::sendToMaster<char>(char* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::sendToMaster<int>(int* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::sendToMaster<float>(float* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::sendToMaster<double>(double* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::iRecv<bool>(bool *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_BYTE, source, tag, comm, request);
  }
}

template <>
void MpiManager::iRecv<char>(char *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, comm, request);
  }
}

template <>
void MpiManager::iRecv<int>(int *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_INT, source, tag, comm, request);
  }
}

template <>
void MpiManager::iRecv<float>(float *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, comm, request);
  }
}

template <>
void MpiManager::iRecv<double>(double *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, comm, request);
  }
}

template <>
void MpiManager::sendRecv<bool>
(bool *sendBuf, bool *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_BYTE, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_BYTE, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<char>
(char *sendBuf, char *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_CHAR, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_CHAR, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<int>
(int *sendBuf, int *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_INT, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_INT, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<float>
(float *sendBuf, float *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_FLOAT, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_FLOAT, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<double>
(double *sendBuf, double *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_DOUBLE, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_DOUBLE, source, tag, comm, &status);
}

template <>
void MpiManager::scatterv_impl<bool>(bool* sendBuf, int* sendCounts, int* displs,
                                     bool* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_BYTE,
               static_cast<void*>(recvBuf),
               recvCount, MPI_BYTE, root, comm);
}

template <>
void MpiManager::scatterv_impl<char>(char* sendBuf, int* sendCounts, int* displs,
                                     char* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_CHAR,
               static_cast<void*>(recvBuf),
               recvCount, MPI_CHAR, root, comm);
}

template <>
void MpiManager::scatterv_impl<int>(int *sendBuf, int* sendCounts, int* displs,
                                    int* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_INT,
               static_cast<void*>(recvBuf),
               recvCount, MPI_INT, root, comm);
}

template <>
void MpiManager::scatterv_impl<float>(float *sendBuf, int* sendCounts, int* displs,
                                      float* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_FLOAT,
               static_cast<void*>(recvBuf),
               recvCount, MPI_FLOAT, root, comm);
}

template <>
void MpiManager::scatterv_impl<double>(double *sendBuf, int* sendCounts, int* displs,
                                       double* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_DOUBLE,
               static_cast<void*>(recvBuf),
               recvCount, MPI_DOUBLE, root, comm);
}

template <>
void MpiManager::gatherv_impl<bool>(bool* sendBuf, int sendCount,
                                    bool* recvBuf, int* recvCounts, int* displs,
                                    int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_BYTE,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_BYTE,
              root, comm);
}

template <>
void MpiManager::gatherv_impl<char>(char* sendBuf, int sendCount,
                                    char* recvBuf, int* recvCounts, int* displs,
                                    int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_CHAR,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_CHAR,
              root, comm);
}

template <>
void MpiManager::gatherv_impl<int>(int* sendBuf, int sendCount,
                                   int* recvBuf, int* recvCounts, int* displs,
                                   int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_INT,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_INT,
              root, comm);
}

template <>
void MpiManager::gatherv_impl<float>(float* sendBuf, int sendCount,
                                     float* recvBuf, int* recvCounts, int* displs,
                                     int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_FLOAT,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_FLOAT,
              root, comm);
}

template <>
void MpiManager::gatherv_impl<double>(double* sendBuf, int sendCount,
                                      double* recvBuf, int* recvCounts, int* displs,
                                      int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_DOUBLE,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_DOUBLE,
              root, comm);
}

template <>
void MpiManager::bCast<bool>(bool* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_BYTE, root, comm);
}

template <>
void MpiManager::bCast<char>(char* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_CHAR, root, comm);
}

template <>
void MpiManager::bCast<int>(int* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_INT, root, comm);
}

template <>
void MpiManager::bCast<float>(float* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_FLOAT, root, comm);
}

template <>
void MpiManager::bCast<double>(double* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_DOUBLE, root, comm);
}

template <>
void MpiManager::bCast<std::string>(std::string* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  int length = (int) sendBuf->size();
  MPI_Bcast(static_cast<void*>(&length), 1, MPI_INT, root, comm);
  char* buffer = new char[length+1];
  if (getRank()==root) {
    std::copy(sendBuf->c_str(), sendBuf->c_str()+length+1, buffer);
  }
  MPI_Bcast(static_cast<void*>(buffer), length+1, MPI_CHAR, root, comm);
  if (getRank()!=root) {
    *sendBuf = buffer;
  }
  delete [] buffer;
}

void MpiManager::bCast(BlockData2D<double,double>& sendData, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendData.getRawData()),
            sendData.getDataSize(), MPI_DOUBLE, root, comm);
}

template <>
void MpiManager::bCastThroughMaster<bool>(bool* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<char>(char* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<int>(int* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<float>(float* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<double>(double* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::reduce<bool>(bool& sendVal, bool& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_BYTE, op, root, comm);
}

template <>
void MpiManager::reduce<char>(char& sendVal, char& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_CHAR, op, root, comm);
}

template <>
void MpiManager::reduce<int>(int& sendVal, int& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_INT, op, root, comm);
}

template <>
void MpiManager::reduce<float>(float& sendVal, float& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_FLOAT, op, root, comm);
}

template <>
void MpiManager::reduce<double>(double& sendVal, double& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_DOUBLE, op, root, comm);
}


template <>
void MpiManager::reduce<BlockData2D<double,int> >(BlockData2D<double,int>& sendVal, BlockData2D<double,int>& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(sendVal.getRawData()),
             static_cast<void*>(recvVal.getRawData()),
             sendVal.getDataSize(), MPI_DOUBLE, op, root, comm);
}

template <>
void MpiManager::reduce<BlockData2D<double,double> >(BlockData2D<double,double>& sendVal, BlockData2D<double,double>& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(sendVal.getRawData()),
             static_cast<void*>(recvVal.getRawData()),
             sendVal.getDataSize(), MPI_DOUBLE, op, root, comm);
}

/*template <>
void MpiManager::reduceVect<bool>(std::vector<bool>& sendVal, std::vector<bool>& recvVal,
                                  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) return;
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_BYTE, op, root, comm);
}
*/
template <>
void MpiManager::reduceVect<char>(std::vector<char>& sendVal, std::vector<char>& recvVal,
                                  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_CHAR, op, root, comm);
}

template <>
void MpiManager::reduceVect<int>(std::vector<int>& sendVal, std::vector<int>& recvVal,
                                 MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_INT, op, root, comm);
}

template <>
void MpiManager::reduceVect<float>(std::vector<float>& sendVal, std::vector<float>& recvVal,
                                   MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_FLOAT, op, root, comm);
}

template <>
void MpiManager::reduceVect<double>(std::vector<double>& sendVal, std::vector<double>& recvVal,
                                    MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_DOUBLE, op, root, comm);
}

template <>
void MpiManager::reduceAndBcast<bool>(bool& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  char recvVal;
  MPI_Reduce(static_cast<void*>(&reductVal), static_cast<void*>(&recvVal), 1, MPI_BYTE, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_BYTE, root, comm);

}

template <>
void MpiManager::reduceAndBcast<char>(char& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  char recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_CHAR, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_CHAR, root, comm);

}

template <>
void MpiManager::reduceAndBcast<int>(int& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  int recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_INT, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_INT, root, comm);

}

template <>
void MpiManager::reduceAndBcast<float>(float& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  float recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_FLOAT, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_FLOAT, root, comm);

}

template <>
void MpiManager::reduceAndBcast<double>(double& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  double recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_DOUBLE, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_DOUBLE, root, comm);

}

template <>
void MpiManager::reduceAndBcast<long>(long& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  long recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_LONG, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_LONG, root, comm);

}

template <>
void MpiManager::reduceAndBcast<unsigned long>(unsigned long& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  unsigned long recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_UNSIGNED_LONG, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_UNSIGNED_LONG, root, comm);

}

template <>
void MpiManager::reduceAndBcast<BlockData2D<double,double>>(BlockData2D<double,double>& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }

  BlockData2D<double,double> recvVal(reductVal.getNx(), reductVal.getNy());

  MPI_Reduce(static_cast<void*>(reductVal.getRawData()),
             static_cast<void*>(recvVal.getRawData()),
             reductVal.getDataSize(), MPI_DOUBLE, op, root, comm);

  reductVal.swap(recvVal);

  MPI_Bcast(static_cast<void*>(reductVal.getRawData()),
            reductVal.getDataSize(), MPI_DOUBLE, root, comm);
}

void MpiManager::wait(MPI_Request* request, MPI_Status* status)
{
  if (!ok) {
    return;
  }
  MPI_Wait(request, status);
}

void MpiManager::waitAll(MpiNonBlockingHelper& mpiNbHelper)
{
  if (!ok) {
    return;
  }
  MPI_Waitall(mpiNbHelper.get_size(), mpiNbHelper.get_mpiRequest(), mpiNbHelper.get_mpiStatus());
}


MpiNonBlockingHelper::MpiNonBlockingHelper()
{
  _size = 0;
}

MpiNonBlockingHelper::MpiNonBlockingHelper(
  MpiNonBlockingHelper const& rhs )
{
  _size          = rhs._size;
  if (_size!=0) {
    allocate(_size);
    for (unsigned i=0; i<_size; i++) {
      _mpiRequest[i] = rhs._mpiRequest[i];
      _mpiStatus[i]  = rhs._mpiStatus[i];
    }
  }
}

MpiNonBlockingHelper MpiNonBlockingHelper::operator= (
  MpiNonBlockingHelper rhs )
{
  MpiNonBlockingHelper tmp(rhs);
  return tmp;
}

void MpiNonBlockingHelper::swap ( MpiNonBlockingHelper& rhs )
{
  std::swap(_size, rhs._size);
  std::swap(_mpiRequest, rhs._mpiRequest);
  std::swap(_mpiStatus, rhs._mpiStatus);
}

MpiNonBlockingHelper::~MpiNonBlockingHelper()
{
  free();
}

void MpiNonBlockingHelper::allocate(unsigned i)
{
  free();
  _size = i;
  _mpiRequest = new MPI_Request [i];
  _mpiStatus  = new MPI_Status [i];
}

void MpiNonBlockingHelper::free()
{
  if (_size!=0) {
    delete [] _mpiRequest;
    delete [] _mpiStatus;
    _size = 0;
  }
}

unsigned const& MpiNonBlockingHelper::get_size() const
{
  return _size;
}

MPI_Request* MpiNonBlockingHelper::get_mpiRequest() const
{
  return _mpiRequest;
}

MPI_Request* MpiNonBlockingHelper::get_mpiRequest(int i) const
{
  assert(size_t(i) < _size);
  return &_mpiRequest[i];
}

MPI_Status* MpiNonBlockingHelper::get_mpiStatus() const
{
  return _mpiStatus;
}


}  // namespace singleton


}

#endif
