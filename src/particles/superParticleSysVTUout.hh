/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn
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

#ifndef SUPERPARTICLESYSVTUOUT_HH
#define SUPERPARTICLESYSVTUOUT_HH

#include "superParticleSysVTUout.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSysVtuWriter<T, PARTICLETYPE>::numofpsys()
{
  return _psys._pSystems.size();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::set(unsigned short pref)
{
  _properties |= pref;
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSysVtuWriter<T, PARTICLETYPE>::SuperParticleSysVtuWriter(const SuperParticleSysVtuWriter<T, PARTICLETYPE>& rhs)
  : _psys(rhs._psys), _name(rhs._name), _properties(rhs._properties), _binary(rhs._binary), _haveMaster(rhs._haveMaster), clout(std::cout, "SuperParticleSysVtuWriter")
{
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSysVtuWriter<T, PARTICLETYPE>::SuperParticleSysVtuWriter(const SuperParticleSysVtuWriter<T, PARTICLETYPE>&& rhs)
  : _psys(rhs._psys), _name(rhs._name), _properties(rhs._properties), _binary(rhs._binary), _haveMaster(rhs._haveMaster), clout(std::cout, "SuperParticleSysVtuWriter")
{
}


template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSysVtuWriter<T, PARTICLETYPE>::SuperParticleSysVtuWriter( SuperParticleSystem3D<T, PARTICLETYPE>& psys,
    std::string const filename, unsigned short properties, bool binary)
  : _psys(psys), _name(filename), _properties(properties), _binary(binary), _haveMaster(false), clout(std::cout, "SuperParticleSysVtuWriter")
{
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::write(int iT)
{
  //std::cout << "Write base" << std::endl;
  int rank = 0;
  int size = 1;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
  size = singleton::mpi().getSize();
#endif

  if (rank == 0) { // master only
    if (!_haveMaster) {
      createMasterFile();
    }

    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                    + createFileName(_name) + "_master.pvd";
    std::string fullNamePVD = singleton::directories().getVtkOutDir() + "data/"
                              + createFileName(_name, iT) + ".pvd";

    preamblePVD(fullNamePVD);         // timestep
    for (int iR = 0; iR < size; iR++) { // cuboid
      std::string namePiece =  "data/" + createFileName(_name, iT, iR) + ".vtu";
      // puts name of .vti piece to a .pvd file [fullNamePVD]
      dataPVD(iT, iR, fullNamePVD, namePiece);
      // adds a namePiece to master.pvd file.
      // To do so we overwrite closePVD() and add new entry.
      dataPVDmaster(iT, iR, fullNamePVDmaster, namePiece);
    } // cuboid
    closePVD(fullNamePVD);            // timestep
  } // master only

  std::string fullNameVTU = singleton::directories().getVtkOutDir()
                            + "data/" + createFileName(_name, iT, rank) + ".vtu";
  preambleVTU(fullNameVTU);
  if (_binary) {
    this->dataArrayBinary(fullNameVTU);
  } else {
    this->dataArray(fullNameVTU);
  }
  closeVTU(fullNameVTU);
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::createMasterFile()
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if (rank == 0) {
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                    + createFileName(_name) + "_master.pvd";
    preamblePVD(fullNamePVDmaster);
    closePVD(fullNamePVDmaster);
    _haveMaster = true;
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::preambleVTU(
  const std::string& fullName)
{
  std::ofstream fout(fullName.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }
  fout << "<?xml version=\"1.0\"?>" << std::endl << std::flush;
  fout
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
      << std::endl;
  fout << "<UnstructuredGrid>" << std::endl;
  fout << "<Piece NumberOfPoints=\"" << _psys.rankNumOfParticles()
       << "\" NumberOfCells=\"" << _psys.rankNumOfParticles() << "\">"
       << std::endl;
  fout << "<PointData Vectors=\"Particles\">" << std::endl;
  fout.close();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::closeVTU(
  const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</UnstructuredGrid>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::preamblePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }

  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n" << "<Collection>\n";
  fout.close();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::closePVD(
  const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::dataPVD(int iT, int iC,
    const std::string& fullNamePVD, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }

  fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" " << iC
       << "\" " << "file=\"" << namePiece << "\"/>\n";
  fout.close();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::dataPVDmaster(int iT, int iC,
    const std::string& fullNamePVDMaster, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVDMaster.c_str(),
                     std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25, std::ios::end); // jump -25 form the end of file to overwrite closePVD

    fout << "<DataSet timestep=\"" << iT << "\" " << "group=\"\" part=\" "
         << iC << "\" " << "file=\"" << namePiece << "\"/>\n";
    fout.close();
    closePVD(fullNamePVDMaster);
  } else {
    clout << "Error: could not open " << fullNamePVDMaster << std::endl;
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::dataArray(
  const std::string& fullName)
{
  //std::cout<< "Base member accessed" << std::endl;
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  if (_properties & particleProperties::radius) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Radius\" NumberOfComponents=\"1\" format=\"ascii\">"
        << std::endl;
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getRad() << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties & particleProperties::mass) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Mass\" NumberOfComponents=\"1\" format=\"ascii\">"
        << std::endl;
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getMass() << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties & particleProperties::cuboid) {
    fout
        << "<DataArray type=\"Int16\" Name=\"Cuboid\" NumberOfComponents=\"1\" format=\"ascii\">"
        << std::endl;
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getCuboid() << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties & particleProperties::active) {
    fout
        << "<DataArray type=\"Int16\" Name=\"Active\" NumberOfComponents=\"1\" format=\"ascii\">"
        << std::endl;
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        if (p.getActive()) {
          fout << "1 ";
        } else {
          fout << "0 ";
        }
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties & particleProperties::velocity) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">"
        << std::endl;
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getVel()[0] << " " << p.getVel()[1] << " " << p.getVel()[2]
             << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties & particleProperties::force) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Force\" NumberOfComponents=\"3\" format=\"ascii\">"
        << std::endl;
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getForce()[0] << " " << p.getForce()[1] << " " << p.getForce()[2]
             << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  fout << "</PointData>" << std::endl;

  fout << "<CellData /> " << std::endl;
  fout << "<Cells>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
       << std::endl;
  int32_t i = 0;
  for (auto pS : _psys._pSystems) {
    for (unsigned int p=0; p<pS->_particles.size(); p++) {
      fout << i++ << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
       << std::endl;
  i = 1;
  for (auto pS : _psys._pSystems) {
    for (unsigned int p=0; p<pS->_particles.size(); p++) {
      fout << i++ << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
       << std::endl;
  for (auto pS : _psys._pSystems) {
    for (unsigned int p=0; p<pS->_particles.size(); p++) {
      fout << 1 << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "</Cells>" << std::endl;
  fout << "<Points>" << std::endl;
  fout
      << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\">"
      << std::endl;

  for (auto pS : _psys._pSystems) {
    for (auto& p : pS->_particles) {
      fout << p.getPos()[0] << " " << p.getPos()[1] << " " << p.getPos()[2] << " ";
    }
  }

  fout << "</DataArray>" << std::endl;
  fout << "</Points>" << std::endl;
  fout << "</Piece>" << std::endl;
  fout.close();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSysVtuWriter<T, PARTICLETYPE>::dataArrayBinary(
  const std::string& fullName)
{
  //std::cout<< "Base member accessed - binary" << std::endl;
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  if (_properties & particleProperties::radius) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Radius\" NumberOfComponents=\"1\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = _psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr, fullSize);
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        const float tmp = float(p.getRad());
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();

    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }

  if (_properties & particleProperties::mass) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Mass\" NumberOfComponents=\"1\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = _psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr, fullSize);
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        const float tmp = float(p.getMass());
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();

    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  if (_properties & particleProperties::cuboid) {
    fout
        << "<DataArray type=\"Int32\" Name=\"Cuboid\" NumberOfComponents=\"1\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = _psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<int> dataEncoder(ofstr,
                                   fullSize);
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        const int tmp = int(p.getCuboid());
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();

    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  if (_properties & particleProperties::active) {
    fout
        << "<DataArray type=\"Int32\" Name=\"Active\" NumberOfComponents=\"1\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = _psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<int> dataEncoder(ofstr, fullSize);
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        int tmp = 0;
        if (p.getActive()) {
          tmp = 1;
        }
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();

    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }

  if (_properties & particleProperties::velocity) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = _psys.rankNumOfParticles() * 3; //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr, fullSize);
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        for (int iDim = 0; iDim < 3; ++iDim) {
          const float tmp = float(p.getVel()[iDim]);
          dataEncoder.encode(&tmp, 1);
        }
      }
    }
    ofstr.close();
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  if (_properties & particleProperties::force) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Force\" NumberOfComponents=\"3\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = _psys.rankNumOfParticles() * 3; //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr, fullSize);
    for (auto pS : _psys._pSystems) {
      for (auto& p : pS->_particles) {
        for (int iDim = 0; iDim < 3; ++iDim) {
          const float tmp = float(p.getForce()[iDim]);
          dataEncoder.encode(&tmp, 1);
        }
      }
    }
    ofstr.close();
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  fout << "</PointData>" << std::endl;

  fout << "<CellData /> " << std::endl;
  fout << "<Cells>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\" encoding=\"base64\">" << std::endl;
  fout.close();
  {
    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = _psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    Base64Encoder<int32_t> dataEncoder(ofstr, fullSize);
    int i = 0;
    for (auto pS : _psys._pSystems) {
      for (unsigned int p=0; p<pS->_particles.size(); p++) {
        const int32_t tmp = i++;
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();
  }
  fout.open(fullName.c_str(), std::ios::out | std::ios::app);
  fout << "</DataArray>" << std::endl;

  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\" encoding=\"base64\">"
       << std::endl;
  fout.close();
  {
    std::ofstream ofstr(fullName.c_str(), std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = _psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    Base64Encoder<int32_t> dataEncoder(ofstr, fullSize);
    int i = 1;
    for (auto pS : _psys._pSystems) {
      for (unsigned int p=0; p<pS->_particles.size(); p++) {
        const int32_t tmp = i++;
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();
  }
  fout.open(fullName.c_str(), std::ios::out | std::ios::app);
  fout << "</DataArray>" << std::endl;


  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"binary\" encoding=\"base64\">"
       << std::endl;
  fout.close();
  {
    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = _psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    Base64Encoder<uint8_t> dataEncoder(ofstr, fullSize);
    for (auto pS : _psys._pSystems) {
      for (unsigned int p=0; p<pS->_particles.size(); p++) {
        const uint8_t tmp = 1;
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();
  }
  fout.open(fullName.c_str(), std::ios::out | std::ios::app);
  fout << "</DataArray>" << std::endl;
  fout << "</Cells>" << std::endl;
  fout << "<Points>" << std::endl;
  fout << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"binary\" encoding=\"base64\">"
       << std::endl;

  fout.close();

  std::ofstream ofstr(fullName.c_str(),
                      std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
  if (!ofstr) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  size_t fullSize = _psys.rankNumOfParticles() * 3; //  how many numbers to write
  size_t binarySize = size_t(fullSize * sizeof(float));
  // writes first number, which have to be the size(byte) of the following data
  Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
  unsigned int uintBinarySize = (unsigned int) binarySize;
  sizeEncoder.encode(&uintBinarySize, 1);
  //  write numbers from functor
  Base64Encoder<float> dataEncoder(ofstr, fullSize);
  for (auto pS : _psys._pSystems) {
    for (auto& p : pS->_particles) {
      for (int iDim = 0; iDim < 3; ++iDim) {
        const float tmp = float(p.getPos()[iDim]);
        dataEncoder.encode(&tmp, 1);
      }
    }
  }
  ofstr.close();
  fout.open(fullName.c_str(), std::ios::out | std::ios::app);

  fout << "</DataArray>" << std::endl;
  fout << "</Points>" << std::endl;
  fout << "</Piece>" << std::endl;
  fout.close();
}



// specialization for magnetic particle

template<typename T>
SuperParticleSysVtuWriterMag<T>::SuperParticleSysVtuWriterMag(const SuperParticleSysVtuWriterMag<T>& rhs) :
  SuperParticleSysVtuWriter<T, MagneticParticle3D>(rhs), _properties(rhs._properties) {}

template<typename T>
SuperParticleSysVtuWriterMag<T>::SuperParticleSysVtuWriterMag(const SuperParticleSysVtuWriterMag<T>&& rhs) :
  SuperParticleSysVtuWriter<T, MagneticParticle3D>(rhs), _properties(rhs._properties) {}


template<typename T>
SuperParticleSysVtuWriterMag<T>::SuperParticleSysVtuWriterMag(
  SuperParticleSystem3D<T, MagneticParticle3D>& psys,
  std::string const filename, const std::bitset<9>& properties, bool binary) :
  SuperParticleSysVtuWriter<T, MagneticParticle3D>(psys, filename, 0, binary), _properties(properties) {}


template<typename T>
void SuperParticleSysVtuWriterMag<T>::dataArrayBinary(
  const std::string& fullName)
{
  //std::cout<< "Special member accessed - binary" << std::endl;
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    this->clout << "Error: could not open " << fullName << std::endl;
  }

  if (_properties.test(pPropRadius)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Radius\" NumberOfComponents=\"1\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr, fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        const float tmp = float(p.getRad());
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();

    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }

  if (_properties.test(pPropMass)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Mass\" NumberOfComponents=\"1\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr, fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        const float tmp = float(p.getMass());
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();

    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropCuboid)) {
    fout
        << "<DataArray type=\"Int32\" Name=\"Cuboid\" NumberOfComponents=\"1\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<int> dataEncoder(ofstr, fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        const int tmp = int(p.getCuboid());
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();

    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropActive)) {
    fout
        << "<DataArray type=\"Int32\" Name=\"Active\" NumberOfComponents=\"1\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<int> dataEncoder(ofstr, fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        int tmp = 0;
        if (p.getActive()) {
          tmp = 1;
        }
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();

    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }

  if (_properties.test(pPropVelocity)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles() * 3; //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr, fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        for (int iDim = 0; iDim < 3; ++iDim) {
          const float tmp = float(p.getVel()[iDim]);
          dataEncoder.encode(&tmp, 1);
        }
      }
    }
    ofstr.close();
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }

  if (_properties.test(pPropForce)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Force\" NumberOfComponents=\"3\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles() * 3; //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr, fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        for (int iDim = 0; iDim < 3; ++iDim) {
          const float tmp = float(p.getForce()[iDim]);
          dataEncoder.encode(&tmp, 1);
        }
      }
    }
    ofstr.close();
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropMoment)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Moment\" NumberOfComponents=\"3\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles() * 3; //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr, fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        for (int iDim = 0; iDim < 3; ++iDim) {
          const float tmp = float(p.getMoment()[iDim]);
          dataEncoder.encode(&tmp, 1);
        }
      }
    }
    ofstr.close();
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropAVel)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"AngularVelocity\" NumberOfComponents=\"3\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary);
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles() * 3; //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr,
                                     fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (auto p : pS->_particles) {
        for (int iDim = 0; iDim < 3; ++iDim) {
          const float tmp = float(p.getAVel()[iDim]);
          dataEncoder.encode(&tmp, 1);
        }
      }
    }
    ofstr.close();
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropTorque)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Torque\" NumberOfComponents=\"3\" format=\"binary\" encoding=\"base64\">"
        << std::endl;
    fout.close();

    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary);
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles() * 3; //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(float));
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float> dataEncoder(ofstr,
                                     fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (auto p : pS->_particles) {
        for (int iDim = 0; iDim < 3; ++iDim) {
          const float tmp = float(p.getTorque()[iDim]);
          dataEncoder.encode(&tmp, 1);
        }
      }
    }
    ofstr.close();
    fout.open(fullName.c_str(), std::ios::out | std::ios::app);
    fout << "</DataArray>" << std::endl;
  }
  fout << "</PointData>" << std::endl;

  fout << "<CellData /> " << std::endl;
  fout << "<Cells>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\" encoding=\"base64\">" << std::endl;
  fout.close();
  {
    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    Base64Encoder<int32_t> dataEncoder(ofstr,
                                       fullSize);
    int i = 0;
    for (auto pS : this->_psys._pSystems) {
      for (unsigned int p=0; p<pS->_particles.size(); p++) {
        const int32_t tmp = i++;
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();
  }
  fout.open(fullName.c_str(), std::ios::out | std::ios::app);
  fout << "</DataArray>" << std::endl;

  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\" encoding=\"base64\">"
       << std::endl;
  fout.close();
  {
    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->_psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    Base64Encoder<int32_t> dataEncoder(ofstr, fullSize);
    int i = 1;
    for (auto pS : this->_psys._pSystems) {
      for (unsigned int p=0; p<pS->_particles.size(); p++) {
        const int32_t tmp = i++;
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();
  }
  fout.open(fullName.c_str(), std::ios::out | std::ios::app);
  fout << "</DataArray>" << std::endl;


  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"binary\" encoding=\"base64\">"
       << std::endl;
  fout.close();
  {
    std::ofstream ofstr(fullName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
    if (!ofstr) {
      this->clout << "Error: could not open " << fullName << std::endl;
    }

    size_t fullSize = this->  _psys.rankNumOfParticles(); //  how many numbers to write
    size_t binarySize = size_t(fullSize * sizeof(int));
    Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
    unsigned int uintBinarySize = (unsigned int) binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    Base64Encoder<uint8_t> dataEncoder(ofstr, fullSize);
    for (auto pS : this->_psys._pSystems) {
      for (unsigned int p=0; p<pS->_particles.size(); p++) {
        const uint8_t tmp = 1;
        dataEncoder.encode(&tmp, 1);
      }
    }
    ofstr.close();
  }
  fout.open(fullName.c_str(), std::ios::out | std::ios::app);
  fout << "</DataArray>" << std::endl;
  fout << "</Cells>" << std::endl;
  fout << "<Points>" << std::endl;
  fout << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"binary\" encoding=\"base64\">"
       << std::endl;

  fout.close();

  std::ofstream ofstr(fullName.c_str(),
                      std::ios::out | std::ios::app | std::ios::binary); // only used for binary output // passed to Base64Encoder
  if (!ofstr) {
    this->clout << "Error: could not open " << fullName << std::endl;
  }

  size_t fullSize = this->_psys.rankNumOfParticles() * 3; //  how many numbers to write
  size_t binarySize = size_t(fullSize * sizeof(float));
  // writes first number, which have to be the size(byte) of the following data
  Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
  unsigned int uintBinarySize = (unsigned int) binarySize;
  sizeEncoder.encode(&uintBinarySize, 1);
  //  write numbers from functor
  Base64Encoder<float> dataEncoder(ofstr, fullSize);
  for (auto pS : this->_psys._pSystems) {
    for (auto& p : pS->_particles) {
      for (int iDim = 0; iDim < 3; ++iDim) {
        const float tmp = float(p.getPos()[iDim]);
        dataEncoder.encode(&tmp, 1);
      }
    }
  }
  ofstr.close();
  fout.open(fullName.c_str(), std::ios::out | std::ios::app);

  fout << "</DataArray>" << std::endl;
  fout << "</Points>" << std::endl;
  fout << "</Piece>" << std::endl;
  fout.close();
}

template<typename T>
void SuperParticleSysVtuWriterMag<T>::dataArray(
  const std::string& fullName)
{
  std::cout<< "Special member accessed" << std::endl;
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    this->clout << "Error: could not open " << fullName << std::endl;
  }

  if (_properties.test(pPropRadius)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Radius\" NumberOfComponents=\"1\" format=\"ascii\">"
        << std::endl;
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getRad() << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropMass)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Mass\" NumberOfComponents=\"1\" format=\"ascii\">"
        << std::endl;
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getMass() << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropCuboid)) {
    fout
        << "<DataArray type=\"Int16\" Name=\"Cuboid\" NumberOfComponents=\"1\" format=\"ascii\">"
        << std::endl;
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getCuboid() << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropActive)) {
    fout
        << "<DataArray type=\"Int16\" Name=\"Active\" NumberOfComponents=\"1\" format=\"ascii\">"
        << std::endl;
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        if (p.getActive()) {
          fout << "1 ";
        } else {
          fout << "0 ";
        }
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropVelocity)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">"
        << std::endl;
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getVel()[0] << " " << p.getVel()[1] << " " << p.getVel()[2]
             << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropForce)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Force\" NumberOfComponents=\"3\" format=\"ascii\">"
        << std::endl;
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getForce()[0] << " " << p.getForce()[1] << " " << p.getForce()[2]
             << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropMoment)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Moment\" NumberOfComponents=\"3\" format=\"ascii\">"
        << std::endl;
    for (auto pS : this->_psys._pSystems) {
      for (auto& p : pS->_particles) {
        fout << p.getMoment()[0] << " " << p.getMoment()[1] << " " << p.getMoment()[2]
             << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  if (_properties.test(pPropTorque)) {
    fout
        << "<DataArray type=\"Float32\" Name=\"Torque\" NumberOfComponents=\"3\" format=\"ascii\">"
        << std::endl;
    for (auto pS : this->_psys._pSystems) {
      for (auto p : pS->_particles) {
        fout << p.getTorque()[0] << " " << p.getTorque()[1] << " " << p.getTorque()[2]
             << " ";
      }
    }
    fout << "</DataArray>" << std::endl;
  }
  fout << "</PointData>" << std::endl;

  fout << "<CellData /> " << std::endl;
  fout << "<Cells>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
       << std::endl;
  int32_t i = 0;
  for (auto pS : this->_psys._pSystems) {
    for (unsigned int p=0; p<pS->_particles.size(); p++) {
      fout << i++ << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
       << std::endl;
  i = 1;
  for (auto pS : this->_psys._pSystems) {
    for (unsigned int p=0; p<pS->_particles.size(); p++) {
      fout << i++ << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
       << std::endl;
  for (auto pS : this->_psys._pSystems) {
    for (unsigned int p=0; p<pS->_particles.size(); p++) {
      fout << 1 << " ";
    }
  }
  fout << "</DataArray>" << std::endl;
  fout << "</Cells>" << std::endl;
  fout << "<Points>" << std::endl;
  fout
      << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\">"
      << std::endl;

  for (auto pS : this->_psys._pSystems) {
    for (auto& p : pS->_particles) {
      fout << p.getPos()[0] << " " << p.getPos()[1] << " " << p.getPos()[2] << " ";
    }
  }

  fout << "</DataArray>" << std::endl;
  fout << "</Points>" << std::endl;
  fout << "</Piece>" << std::endl;
  fout.close();
}

template<typename T>
void SuperParticleSysVtuWriterMag<T>::write(int iT)
{
  //std::cout << "Write derived" << std::endl;
  int rank = 0;
  int size = 1;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
  size = singleton::mpi().getSize();
#endif

  if (rank == 0) { // master only
    if (!this->_haveMaster) {
      this->createMasterFile();
    }

    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                    + createFileName(this->_name) + "_master.pvd";
    std::string fullNamePVD = singleton::directories().getVtkOutDir() + "data/"
                              + createFileName(this->_name, iT) + ".pvd";

    this->preamblePVD(fullNamePVD);         // timestep
    for (int iR = 0; iR < size; iR++) { // cuboid
      std::string namePiece =  "data/" + createFileName(this->_name, iT, iR) + ".vtu";
      // puts name of .vti piece to a .pvd file [fullNamePVD]
      this->dataPVD(iT, iR, fullNamePVD, namePiece);
      // adds a namePiece to master.pvd file.
      // To do so we overwrite closePVD() and add new entry.
      this->dataPVDmaster(iT, iR, fullNamePVDmaster, namePiece);
    } // cuboid
    this->closePVD(fullNamePVD);            // timestep
  } // master only

  std::string fullNameVTU = singleton::directories().getVtkOutDir()
                            + "data/" + createFileName(this->_name, iT, rank) + ".vtu";
  this->preambleVTU(fullNameVTU);
  if (this->_binary) {
    this->dataArrayBinary(fullNameVTU);
  } else {
    this->dataArray(fullNameVTU);
  }
  this->closeVTU(fullNameVTU);
}

template<typename T>
void SuperParticleSysVtuWriterMag<T>::set(int pref)
{
  _properties.set(pref);
}



}  // namespace OLB

#endif /* SUPERPARTICLESYSVTUOUT_HH */
