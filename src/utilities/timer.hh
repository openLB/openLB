/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011 Lukas Baron, Mathias J. Krause
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

#ifndef TIMER_HH
#define TIMER_HH

#include "timer.h"
#include "communication/mpiManager.h"

namespace olb {

namespace util {

template<typename T>
Timer<T>::Timer(int maxTimeSteps, size_t numFluidCells)
  : clout(std::cout,"Timer"),
    deltaTS(0), curTS(0), maxTS(maxTimeSteps),
    numFC(numFluidCells), rtRemMs(1) // avoids some stupid numbers in first call of printStep() (not for T=double)

{
  tp = nullptr;
}

template<typename T>
T Timer<T>::timevalDiffTimeMs(timeval end, timeval start)
{
  T msDiff;
  msDiff = 1000*(end.tv_sec - start.tv_sec)
           +(end.tv_usec-start.tv_usec)/1000;
  return std::max<T>(msDiff, 1.0);
}

template<typename T>
T Timer<T>::getMLUPs()
{
  T mlups = (numFC * deltaTS) / (timevalDiffTimeMs(msTimeCur, msTimeLast)*1000);
  return mlups;
}

template<typename T>
T Timer<T>::getMLUPps()
{
  T mlupps = getMLUPs()/singleton::mpi().getSize();
  return mlupps;
}

template<typename T>
T Timer<T>::getTotalMLUPs()
{
  T tmlups = ((T)numFC * curTS) / (timevalDiffTimeMs(msTimeEnd, msTimeStart)*1000);
  return tmlups;
}

template<typename T>
T Timer<T>::getTotalMLUPps()
{
  T tmlupps = getTotalMLUPs()/singleton::mpi().getSize();
  return tmlupps;
}


template<typename T>
void Timer<T>::start()
{
  sTimeStart = time(tp);          // time in s
  gettimeofday(&msTimeStart, nullptr);  // time in ms
  gettimeofday(&msTimeCur, nullptr);    // time in ms, here only necessary for MLUP-calculations
  cpuTimeStart = clock();         //cpu-time
}

template<typename T>
void Timer<T>::update(int currentTimeStep)    // Is int sufficient? Is it possible/desirable to have non-integer time steps?
{

  cpuTimeCur = clock();           // CPU-time
  sTimeCur   = time(tp);          // time in s
  msTimeLast = msTimeCur;
  gettimeofday(&msTimeCur, nullptr);    // time in ms

  // calculate and update missing time-values
  deltaTS = currentTimeStep - curTS;      // this makes multiple calls
  curTS   = currentTimeStep;              // of update() critical

  rtPas   = difftime(sTimeCur,sTimeStart);                // here calculation is based on s-time
  rtTot   = rtPas*maxTS/std::max<int>(curTS, 1);
  rtRem   = rtTot-rtPas;

  rtPasMs = timevalDiffTimeMs(msTimeCur, msTimeStart);    // here with ms-time as timeval-value
  rtTotMs = rtPasMs*maxTS/std::max<int>(curTS, 1);
  rtRemMs = rtTotMs-rtPasMs;

  ctPas   = (cpuTimeCur-cpuTimeStart)/CLOCKS_PER_SEC;     // and here the same for CPU-time
  ctTot   = ctPas*maxTS/std::max<int>(curTS, 1);
  ctRem   = ctTot-ctPas;

}

template<typename T>
void Timer<T>::stop()
{
  cpuTimeEnd = clock();           // cpu-time
  sTimeEnd = time(tp);            // time in s
  gettimeofday(&msTimeEnd, nullptr);    // time in ms
}

template<typename T>
double Timer<T>::getTotalCpuTime()
{
  return (cpuTimeEnd-cpuTimeStart)/CLOCKS_PER_SEC;
}

template<typename T>
T Timer<T>::getTotalRealTime()
{
  return difftime(sTimeEnd,sTimeStart);
}

template<typename T>
T Timer<T>::getTotalRealTimeMs()
{
  return timevalDiffTimeMs(msTimeEnd, msTimeStart);
}

template<typename T>
void Timer<T>::print(int currentTimeStep,  int printMode)
{
  if (currentTimeStep!=curTS) {
    update(currentTimeStep);
  }
  printStep(printMode);
}

template<typename T>
void Timer<T>::printStep(int printMode)
{
  switch (printMode) {
  case 0: //single-line layout, usable for data extraction as csv
    clout
        << "step=" << curTS << "; "
        //      << "stepMax=" << maxTS << "; "
        << "percent=" << 100.0*curTS/maxTS << "; "
        << "passedTime=" << (double)rtPasMs/1000 << "; "
        //      << "totalTime=" << (double)rtTotMs/1000 << "; "
        << "remTime=" << rtRemMs/1000 << "; "
        << "MLUPs=" << getMLUPs()
        << std::endl;
    break;

  case 1: //single-line layout (not conform with output-rules)
    clout
        << "latticeTS: "
        << curTS << "/" << maxTS << " (" << 100*curTS/maxTS << "%); "
        << "pas/totTime: "
        << std::setprecision(2) << std::fixed << (double)rtPasMs/1000 << "/"
        << std::setprecision(1) << std::fixed << (double)rtTotMs/1000 << "s; "
        << "remTime: "
        << std::setw(2) << (int)((double)rtRemMs/1000)/60 << "m " << std::setfill('0') << std::setw(4) << (double)((int)((double)rtRemMs/100)%600)/10 << "s; "
        << std::setfill(' ')
        << "MLUPs: " << getMLUPs()
        << std::endl;
    break;

  case 2: //pretty double line layout in colums, but non-conform
    clout
        << std::setw(21) << std::left << "Lattice-Timesteps"
        << std::setw(17) << std::left << "| CPU time/estim"
        << std::setw(18) << std::left << "| REAL time/estim"
        << std::setw(6)  << std::left << "| ETA"
        << std::setw(6)  << std::left << "| MLUPs"
        << std::endl << std::right
        << std::setw(6) << std::setprecision(2) << std::fixed << curTS << "/" << std::setw(6) << maxTS << " (" << std::setw(3) << 100*curTS/maxTS << "%) |"
        << std::setw(7) << ctPas << "/" << std::setw(7) << ctTot << " |"
        << std::setw(8) << (double)rtPasMs/1000 << "/" << std::setw(7) << (double)rtTotMs/1000 << " |"
        << std::setw(4) << (int)rtRemMs/1000+1 << " |"
        << std::setw(6) << getMLUPs()
        << std::endl;
    break;

  case 3: //performance output only
    clout
        << "step " << curTS << "; "
        << "MLUPs=" << std::setw(8) << getMLUPs() << ", MLUPps=" << std::setw(8) << getMLUPps() << std::endl;
    break;

  default:
    clout << "Error in function printStep in class_timer.h: printMode="<<printMode<<" not found" << std::endl << std::flush;
  }
}

template<typename T>
void Timer<T>::printSummary()
{
  clout << std::endl;
  clout << "----------------Summary:Timer----------------" << std::endl;
  clout << "measured time (rt) : " << (int)getTotalRealTimeMs()/1000 << "." << (int)getTotalRealTimeMs()-(int)getTotalRealTimeMs()/1000*1000 << "s" << std::endl;
  clout << "measured time (cpu): " << std::setprecision(3) << std::fixed << getTotalCpuTime() << "s" << std::endl;
  clout << "average MLUPs :       " << getTotalMLUPs()  << std::endl;
  clout << "average MLUPps:       " << getTotalMLUPps() << std::endl;
  clout << "---------------------------------------------" << std::endl;
}

template<typename T>
void Timer<T>::printShortSummary()
{
  clout << "realTime=" << (int)getTotalRealTimeMs()/1000 << "." << (int)getTotalRealTimeMs()-(int)getTotalRealTimeMs()/1000*1000
        << "; cpuTime=" <<  std::setprecision(3) << std::fixed << getTotalCpuTime() << std::endl;
}

// Factory function /////////////////////////////////

template<typename T, typename DESCRIPTOR>
Timer<T>* createTimer(XMLreader& param, const UnitConverter<T,DESCRIPTOR>& converter, size_t numLatticePoints)
{
  OstreamManager clout(std::cout,"createTimer");

  // initialize parameters with some default values
  T physMaxT = T();
  T physStartT = T();

  // fetch xml Data and error handling
  if ( ! param["Application"]["PhysParameters"]["PhysMaxTime"].read(physMaxT) ) {
    if ( ! param["Application"]["PhysParam"]["MaxTime"].read(physStartT) ) {
      clout << "PhysMaxTime not found" << std::endl;
    }
    else
    {
      clout << "Application::PhysParam::MaxTime needs to be renamed to Application::PhysParameters::PhysMaxTime" << std::endl;
    }
  }
//  if ( ! param["Application"]["PhysParam"]["MaxStartTime"].read(physStartT) ) {
//    clout << "PhysStartTime not found" << std::endl;
//  }

  // variable processing according to the constructor
  int maxT = converter.getLatticeTime(physMaxT) + converter.getLatticeTime(physStartT);

  //return some default values that produce reasonable output (e.g.
  // zero); in best case there should be no output at all (TODO)
  return new Timer<T>(maxT, numLatticePoints);
}

} // namespace util

} // namespace olb

#endif
