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

#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include "io/ostreamManager.h"
#include "io/xmlReader.h"
#include "core/unitConverter.h"

namespace olb {

namespace util {

// using namespace olb; //necessary? Problems have occured with singleton::xxx

/** \file
 This class allows calculation and display of various time data including remaining runtime data in cpu and in user time. Output is possible during computation and as a summary after computation. Time can be measured independently from in-computation-calls which allows time statistics without touching the overall computation time of the LBM-algorithm. Moreover it is possible to determine a performance number, the Mega DESCRIPTOR Updates per second (MLUPs).
The template type T denotes the internal representation of time differences (in ms, s and for cpu-time). Reasonable values are eg. double and long.
 Some thoughts about data range: Actually, only positive numbers should occur (=> unsigned). The chosen data format should be capable of representing at least 7 days ( = 604.800.000 ms) of computation time. The biggest signed long32 number ist 2.147.483.647 (~3 1/2 weeks), great! For everything above, use double with it's floating point arithmetic.
*/

/** How to use in Code:
  <pre>
  Timer timer;
  timer.inizialize();

  timer.start();

  loop(i=1..iMax)
  {
    do_some_calculation();
    timer.print(i);
  } end_loop;

  timer.stop();

  timer.printSummary();
  </pre>
*/

/* BUGS:
- cavity2d with 20 seconds yields average MLUPs values much smaller than the ones 'live in computation' (with T=int, due to integer overflow)
*/

/// class for measurement of computation time
template<typename T>
class Timer {

private:
  mutable OstreamManager clout;

  /* INFORMATION about time data types:
      timeval: Struct with "long tv_usec" and "unsigned long tv_sec".
      time_t : seconds since 1st january 1970, 00:00:00 (GMT).
  */

  // parameter for time measurement
  double        cpuTimeStart, cpuTimeCur, cpuTimeEnd;         // in cpu-time
  time_t          sTimeStart,   sTimeCur,   sTimeEnd, *tp;    // in seconds
  timeval        msTimeStart,  msTimeCur,  msTimeEnd;         // in ms
  timeval        msTimeLast;                                  // in ms, for MegaLatticeUpdate-measurements only
  T deltaTS;                                                  // time-step difference since last call of update()

  // Input-parameter of the algorithm to compute remaining runtimes and performance (actually constants)
  int curTS;      // current lattice time step
  int maxTS;      // total number of lattice time steps that are intended to be computed
  size_t numFC;      // number of fluid cells (depending from size and dimension of the domain)

  // parameter-prefixes for output
  /* prefix-explanation:
      rt: realtime-values, elapsed time on a wall clock
      ct: cpu-time values
      lt: lattice-time values, elapsed time within the simulated system
  */

  int    ltTot, ltPas, ltRem;         // lattice time (time steps)
  double ctPas, ctRem, ctTot;         // cpu time
  T      rtPas, rtRem, rtTot;         // times in s
  T      rtPasMs, rtRemMs, rtTotMs;   // times in ms

public:
  /// initializes timer with the given values, abbreviation to Timer() + initialize(int,int)
  Timer(int maxTimeSteps, size_t numFluidCells=1);

  /// returns the time difference between two timeval objects in ms
  /** The timeval data type is used in the variables for ms-time measurement. \sa getTotalRealTimeMs*/
  T timevalDiffTimeMs(timeval end, timeval start);

  /// returns Million DESCRIPTOR Site Updates per second (all processes together)
  T getMLUPs();

  /// returns Million DESCRIPTOR Site Updates per second and process
  T getMLUPps();

  /// returns average Million DESCRIPTOR Site Updates per second between start() and stop()
  T getTotalMLUPs();

  /// returns average Million DESCRIPTOR Site Updates per second and process between start() and stop()
  T getTotalMLUPps();

  /// (Re-)sets start value for time measurement.
  void start();

  /// Updates all time values of interest during computation
  /** i.e. elapsed time, remaining time and total time. Recommended to be used directly before printStep to output up-to-date values.*/
  void update(int currentTimeStep);

  /// Terminates time measurement and sets end value.
  /** It is necessary to call this function immediately after the computation loop has terminated. Otherwise printSummary() would display incorrect values. \sa printSummary()*/
  void stop();

  /// Returns the total cpu time in seconds between start() and stop().
  double getTotalCpuTime();

  /// Returns the total measured time between start() and stop() in seconds.
  T getTotalRealTime();

  /// Returns the total measured time between start() and stop() in ms.
  T getTotalRealTimeMs();

  /// Prints a one-line-summary of the values calculated in update() for use during computation.
  /** Displays timer informations about current calculation with two different output modes. Before calling printStep() one should first recalculate all values with update().
    \param printMode Value for selecting style of output. <br>0 = default semicolon-separated single row mode without cpu-time <br>1 = single row mode without cpu-time <br>2 = nicely formatted two-line layout including current and remaining cpu-time
    \sa update()
  */
  void printStep(int printMode=0);

  /// Performs an update() followed by a printStep().
  /** Automatically calls the function update(currentTimeStep) (if not yet done with currentTimeStep) and displays the timer's during-iteration-information.
    \param currentTimeStep current iteration value, e.g. i or iT<br>
    \param printMode mode of display style passed to printStep()
    \sa printStep()
  */
  void print(int currentTimeStep, int printMode=0);

  /// Prints a (short) summary containing the overall time consumption in real and in cpu time for use after computation.
  void printSummary();
  /// Prints a short summary containing only time consumptions (real and cpu time)
  void printShortSummary();
};

// Factory function /////////////////////////////////
template<typename T, typename DESCRIPTOR>
// Timer<T>* createTimer(XMLreader& param);
Timer<T>* createTimer(XMLreader& param, const UnitConverter<T,DESCRIPTOR>& converter, size_t numLatticePoints);
/////////////////////////////////////////////////////

} // namespace util

} // namespace olb

#endif
