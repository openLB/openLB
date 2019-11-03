/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2007 Jonas Latt,
 *                2015-2019 Mathias J. Krause, Adrian Kummerlaender
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

/** \file
 * Definition of a LB cell -- header file.
 */
#ifndef CELL_H
#define CELL_H

#include "olbDebug.h"
#include "serializer.h"
#include "dynamics/latticeDescriptors.h"
#include "dynamics/dynamics.h"

namespace olb {

template<typename T, typename DESCRIPTORBASE>
class CellBase {
protected:
  /// The lattice populations and fields are stored as a C-array.
  T data[DESCRIPTORBASE::size()]; ///< distribution functions and additional fields

public:
  /// Read-write access to distribution functions.
  /**
   * \param iPop index of the accessed distribution function
   *
   * Note that for legacy-purposes this method allows access to all data of the cell
   * i.e. not just its distribution but also all fields that are declared by the
   * descriptor. Use of e.g. Cell::getFieldPointer or Cell::defineField is strongly
   * recommended when writing new code.
   **/
  T& operator[](int const& iPop)
  {
    OLB_PRECONDITION( iPop < DESCRIPTORBASE::size() );
    return data[iPop];
  }

  /// Read-only access to distribution functions.
  /**
   * \param iPop index of the accessed distribution function
   **/
  T const& operator[](int const& iPop) const
  {
    OLB_PRECONDITION( iPop < DESCRIPTORBASE::size() );
    return data[iPop];
  }
};


/// A LB lattice cell.
/**
 * A cell contains the q values of the distribution functions f on one lattice
 * point, additional "external" fields as well as a pointer to the dynamics of
 * the cell. Thanks to this pointer, one can have a space dependend definition
 * of the dynamics. This mechanism is useful e.g. for the implementation of
 * boundary conditions, or an inhomogeneous body force.
 *
 * The dynamics object is not owned by the class and as such not destructed in
 * the Cell destructor.
 *
 * This class is not intended to be derived from.
 */
template<typename T, typename DESCRIPTOR>
class Cell :
  public CellBase<T, typename DESCRIPTOR::BaseDescriptor>,
  public Serializable {
private:
  Dynamics<T,DESCRIPTOR>* dynamics;  ///< local LB dynamics

public:
  /// Default constructor.
  Cell();
  /// Constructor, to be used whenever possible.
  Cell(Dynamics<T,DESCRIPTOR>* dynamics_);

  /// Return pointer to FIELD of cell
  template <typename FIELD, typename X = DESCRIPTOR>
  utilities::meta::enable_if_t<X::template provides<FIELD>(), T*>
  getFieldPointer()
  {
    const int offset = DESCRIPTOR::template index<FIELD>();
    return &(this->data[offset]);
  }

  template <typename FIELD, typename X = DESCRIPTOR>
  utilities::meta::enable_if_t<!X::template provides<FIELD>(), T*>
  getFieldPointer()
  {
    throw std::invalid_argument("DESCRIPTOR does not provide FIELD.");
    return nullptr;
  }

  /// Return read-only pointer to FIELD of cell
  template <typename FIELD, typename X = DESCRIPTOR>
  utilities::meta::enable_if_t<X::template provides<FIELD>(), const T*>
  getFieldPointer() const
  {
    const int offset = DESCRIPTOR::template index<FIELD>();
    return &(this->data[offset]);
  }

  template <typename FIELD, typename X = DESCRIPTOR>
  utilities::meta::enable_if_t<!X::template provides<FIELD>(), const T*>
  getFieldPointer() const
  {
    throw std::invalid_argument("DESCRIPTOR does not provide FIELD.");
    return nullptr;
  }

  /// Return copy of FIELD as a vector
  template <typename FIELD, typename X = DESCRIPTOR>
  utilities::meta::enable_if_t<(X::template size<FIELD>() > 1), Vector<T,X::template size<FIELD>()>>
  getField() const
  {
    return Vector<T,DESCRIPTOR::template size<FIELD>()>(
      getFieldPointer<FIELD>()
    );
  }

  /// Return copy of FIELD as a scalar
  template <typename FIELD, typename X = DESCRIPTOR>
  utilities::meta::enable_if_t<(X::template size<FIELD>() == 1), T>
  getField() const
  {
    return getFieldPointer<FIELD>()[0];
  }

  /// Set value of FIELD from a vector
  template <typename FIELD, typename X = DESCRIPTOR>
  utilities::meta::enable_if_t<(X::template size<FIELD>() > 1), void>
  setField(const Vector<T,DESCRIPTOR::template size<FIELD>()>& field)
  {
    std::copy_n(
      field.data,
      DESCRIPTOR::template size<FIELD>(),
      getFieldPointer<FIELD>());
  }

  /// Set value of FIELD from a scalar
  template <typename FIELD, typename X = DESCRIPTOR>
  utilities::meta::enable_if_t<(X::template size<FIELD>() == 1), void>
  setField(T value)
  {
    getFieldPointer<FIELD>()[0] = value;
  }

  /// Copy FIELD content to given memory location
  template <typename FIELD>
  void computeField(T* ext) const
  {
    const T* field = getFieldPointer<FIELD>();
    for (int iExt=0; iExt < DESCRIPTOR::template size<FIELD>(); ++iExt) {
      ext[iExt] = field[iExt];
    }
  }

  /// Set FIELD value from given memory location
  template <typename FIELD>
  void defineField(const T* ext)
  {
    T* field = getFieldPointer<FIELD>();
    for (int iExt=0; iExt < DESCRIPTOR::template size<FIELD>(); ++iExt) {
      field[iExt] = ext[iExt];
    }
  }

  /// Add to FIELD from given memory location
  /**
   * Similar to defineField(),but instead of replacing existing values
   * the data at ext is added onto the existing values.
   **/
  template <typename FIELD>
  inline void addField(const T* ext)
  {
    T* field = getFieldPointer<FIELD>();
    for (int iExt=0; iExt < DESCRIPTOR::template size<FIELD>(); ++iExt) {
      field[iExt] += ext[iExt];
    }
  }

  /// Multiply FIELD with values at given memory location
  /**
   * Similar to defineField(), but instead of replacing existing values
   * the data at ext is multiplied to the existing values.
   **/
  template <typename FIELD>
  inline void multiplyField(const T* ext)
  {
    T* field = getFieldPointer<FIELD>();
    for (int iExt=0; iExt < DESCRIPTOR::template size<FIELD>(); ++iExt) {
      field[iExt] *= ext[iExt];
    }
  }

  /// Define or re-define dynamics of the cell.
  /**
   * \param dynamics_ a pointer to the dynamics object, whos memory management
   *                  falls under the responsibility of the user
   **/
  void defineDynamics(Dynamics<T,DESCRIPTOR>* dynamics_);
  /// Get a non-modifiable pointer to the dynamics
  Dynamics<T,DESCRIPTOR> const* getDynamics() const;
  /// Get a non-modifiable pointer to the dynamics
  Dynamics<T,DESCRIPTOR>* getDynamics();

  // The following helper functions forward the function call
  // to the Dynamics object
public:
  /// Apply LB collision to the cell according to local dynamics.
  void collide(LatticeStatistics<T>& statistics)
  {
    OLB_PRECONDITION( dynamics );
    dynamics->collide(*this, statistics);
  }

  /// Compute particle density on the cell.
  /** \return particle density
   */
  T computeRho() const
  {
    OLB_PRECONDITION( dynamics );
    return dynamics->computeRho(*this);
  }
  /// Compute fluid velocity on the cell.
  /** \param u fluid velocity
   */
  void computeU(T u[descriptors::d<DESCRIPTOR>()]) const
  {
    OLB_PRECONDITION( dynamics );
    dynamics->computeU(*this, u);
  }
  /// Compute fluid momentum (j = rho * u) on the cell.
  /** \param j fluid momentum
   */
  void computeJ(T j[descriptors::d<DESCRIPTOR>()]) const
  {
    OLB_PRECONDITION( dynamics );
    dynamics->computeJ(*this, j);
  }
  /// Compute components of the stress tensor on the cell.
  /** \param pi stress tensor */
  void computeStress (
    T pi[util::TensorVal<DESCRIPTOR >::n]) const
  {
    OLB_PRECONDITION( dynamics );
    T rho, u[descriptors::d<DESCRIPTOR>()];
    dynamics->computeRhoU(*this, rho, u);
    dynamics->computeStress(*this, rho, u, pi);
  }
  /// Compute fluid velocity and particle density on the cell.
  /** \param rho particle density
   *  \param u fluid velocity
   */
  void computeRhoU(T& rho, T u[descriptors::d<DESCRIPTOR>()]) const
  {
    OLB_PRECONDITION( dynamics );
    dynamics->computeRhoU(*this, rho, u);
  }
  /// Compute equilibrium part of cell distribution
  void computeFeq(T fEq[descriptors::q<DESCRIPTOR>()]) const;
  /// Compute non-equilibrium part of cell distribution
  void computeFneq(T fNeq[descriptors::q<DESCRIPTOR>()]) const;
  /// Compute all momenta on the celll, up to second order.
  /** \param rho particle density
   *  \param u fluid velocity
   *  \param pi stress tensor
   */
  void computeAllMomenta (
    T& rho, T u[descriptors::d<DESCRIPTOR>()],
    T pi[util::TensorVal<DESCRIPTOR >::n] ) const
  {
    OLB_PRECONDITION( dynamics );
    dynamics->computeAllMomenta(*this, rho, u, pi);
  }
  /// Set particle density on the cell.
  /** \param rho particle density
   */
  void defineRho(T rho)
  {
    OLB_PRECONDITION( dynamics );
    dynamics->defineRho(*this, rho);
  }
  /// Set fluid velocity on the cell.
  /** \param u fluid velocity
   */
  void defineU(const T u[descriptors::d<DESCRIPTOR>()])
  {
    OLB_PRECONDITION( dynamics );
    dynamics->defineU(*this, u);
  }
  /// Define fluid velocity and particle density on the cell.
  /** \param rho particle density
   *  \param u fluid velocity
   */
  void defineRhoU(T rho, const T u[descriptors::d<DESCRIPTOR>()])
  {
    OLB_PRECONDITION( dynamics );
    dynamics->defineRhoU(*this, rho, u);
  }
  /// Define particle populations through the dynamics object.
  /** This method is similar to operator[]: it modifies the
   * value of all the particle populations.
   */
  void definePopulations(const T* f_)
  {
    for (int iPop = 0; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
      this->data[iPop] = f_[iPop];
    }
  }
  /// Initialize all f values to their local equilibrium
  void iniEquilibrium(T rho, const T u[descriptors::d<DESCRIPTOR>()])
  {
    OLB_PRECONDITION( dynamics );
    dynamics->iniEquilibrium(*this, rho, u);
  }
  /// Revert ("bounce-back") the distribution functions.
  void revert();
  void serialize(T* data) const;
  void unSerialize(T const* data);

  /// \return the number of data blocks for the serializable interface
  std::size_t getNblock() const override
  {
    return 1;
  };
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// \return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

private:
  void initializeData();
};

template<typename T, typename DESCRIPTOR>
struct WriteCellFunctional {
  virtual ~WriteCellFunctional() { };
  virtual void apply(Cell<T,DESCRIPTOR>& cell, int pos[descriptors::d<DESCRIPTOR>()]) const =0;
};

}  // namespace olb

#endif
