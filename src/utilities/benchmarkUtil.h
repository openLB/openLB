#ifndef BENCHMARK_UTIL_H
#define BENCHMARK_UTIL_H

#include <deque>
#include "io/ostreamManager.h"
#include "functors/analytical/analyticalF.h"

namespace olb {

namespace util {

/// Check time-convergence of a scalar.
/** This class is useful, for example to check convergence of
 * the velocity field for the simulation of a stationary flow.
 * Convergence is claimed when the standard deviation of the
 * monitored value is smaller than epsilon times the average.
 * The statistics are taken over a macroscopic time scale of the
 * system.
 */
template<typename T>
class ValueTracer {
public:
  /// Ctor.
  /** \param u The characteristic velocity of the system, for
   *          computation of the characteristic time scale.
   * \param L The characteristic length of the system, for
   *          computation of the characteristic time scale.
   * \param epsilon Precision of the convergence.
   */
  ValueTracer(T u, T L, T epsilon);
  /// Ctor.
  /** \param deltaT   corresponds to iteration steps (averaging)
  \*  \param epsilon  allowed derivation of quantity
  */
  ValueTracer(int deltaT, T epsilon);
  /// Change values of u and L to update characteristic scales of the system.
  void resetScale(T u, T L);
  /// reinitializes the values
  void resetValues();
  /// Get characteristic time scale.
  int getDeltaT() const;
  /// Feed the object with a new measured scalar.
  void takeValue(T val, bool doPrint=false);
  /// Test for convergence, with respect to stdDev.
  bool hasConverged() const;
  /// Test for convergence, with respect to difference between min and max value;
  bool convergenceCheck() const;
  /// Test for convergence, with respect to difference between min and max value;
  bool hasConvergedMinMax() const;
  T computeAverage() const;
  T computeStdDev(T average) const;
  void setEpsilon(T epsilon);
private:
  int    _deltaT;
  T      _epsilon;
  int    _t;
  bool   _converged;
  std::deque<T> _values;
  mutable OstreamManager clout;
};

/// Propose successive test values of a scalar (e.g. Re) to check stability of a system.
/** At first, the stability limit is explored by constant
 * increments/decrements of the scalar, and then, by successive
 * bisection.
 */
template<typename T>
class BisectStepper {
public:
  /// The only constructor.
  /** \param _iniVal Initial guess for the stability limit.
   * \param _step   Step size at which the value is initially
   *                incremented/decremented.
   */
  BisectStepper(T _iniVal, T _step=0.);
  /// Get new value, and indicate if the previous value yielded a stable system or not.
  T getVal(bool stable, bool doPrint=false);
  /// Test for convergence.
  bool hasConverged(T epsilon) const;
private:
  T iniVal, currentVal, lowerVal, upperVal;
  T step;
  enum {first, up, down, bisect} state;
  mutable OstreamManager clout;
};

/// 1D Newton simple scheme
template<typename T>
class Newton1D {

protected:
  AnalyticalF1D<T,T>& _f;
  AnalyticalDiffFD1D<T> _df;
  T _yValue;
  T _eps;
  int _maxIterations;

public:
  Newton1D(AnalyticalF1D<T,T>& f, T yValue = T(), T eps = 1.e-8, int maxIterations = 100);

  T solve(T startValue, bool print=false);
};

/// Trapezoidal rule
template<typename T>
class TrapezRuleInt1D {

protected:
  AnalyticalF1D<T,T>& _f;

public:
  TrapezRuleInt1D(AnalyticalF1D<T,T>& f);

  T integrate(T min, T max, int nSteps);
};

/** Simple circular buffer to compute average and other quantities
 *  over pre-defined temporal windows.
 *  Works with every T supporting += (T or scalar), and /= (scalar) operations,
 *  including double and Vector<double, size>.
 */
template<typename T>
class CircularBuffer {
public:
  CircularBuffer(int size);
  /// insert a new entry ed eventually erases the oldest one
  void insert(T entry);
  /// average over all the entries
  T average();
  /// get reference to the last entry for pos=0, the second-to-last for pos=1, and so on
  T& get(int pos);
  /// return size of the buffer
  int getSize();
  
private:
  int _size;
  std::vector<T> _data;
};

} // namespace util

} // namespace olb

#endif
