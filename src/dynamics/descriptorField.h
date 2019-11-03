/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef DESCRIPTOR_FIELD_H
#define DESCRIPTOR_FIELD_H

#include <type_traits>
#include <stdexcept>

namespace olb {

namespace descriptors {

/// \defgroup descriptor
//@{

/// Base of a descriptor field whose size is defined by A*D + B*Q + C
template <unsigned C, unsigned A=0, unsigned B=0>
struct DESCRIPTOR_FIELD_BASE {

  /// Deleted constructor to enforce pure usage as type and prevent implicit narrowing conversions
  DESCRIPTOR_FIELD_BASE() = delete;

  /// Evaluates the size function
  /**
   * To be called by DESCRIPTOR_BASE
   **/
  template <unsigned D, unsigned Q>
  static constexpr unsigned size()
  {
    return A * D + B * Q + C;
  }

  /// Returns global index from local index and provides out_of_range safety
  template <unsigned D, unsigned Q>
  static constexpr unsigned getLocalIndex(const unsigned localIndex)
  {
    return localIndex < (A*D+B*Q+C) ? localIndex : throw std::out_of_range("Index exceeds data field");
  }

};

/// \defgroup descriptor_fields Set of common descriptor fields
/// \ingroup descriptor
//@{

// Field types need to be distinct (i.e. not aliases) in order for `getDescriptorFieldBeginsAt` to work
// (Field size parametrized by: Cs + Ds*D + Qs*Q)          Cs Ds Qs
struct POPULATION           : public DESCRIPTOR_FIELD_BASE<0,  0, 1> { };
struct VELOCITY             : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct VELOCITY2            : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct SOURCE               : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct FORCE                : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct EXTERNAL_FORCE       : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct TAU_EFF              : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct CUTOFF_KIN_ENERGY    : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct CUTOFF_HEAT_FLUX     : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct CHEM_POTENTIAL       : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct V6                   : public DESCRIPTOR_FIELD_BASE<6,  0, 0> { };
struct V12                  : public DESCRIPTOR_FIELD_BASE<12, 0, 0> { };
struct OMEGA                : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct G                    : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct EPSILON              : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct BODY_FORCE           : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct K                    : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct NU                   : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct POROSITY             : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct VELOCITY_NUMERATOR   : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct VELOCITY_DENOMINATOR : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct LOCAL_DRAG           : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct VELOCITY_SOLID       : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct COORDINATE           : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct F                    : public DESCRIPTOR_FIELD_BASE<0,  0, 1> { };
struct DJDF                 : public DESCRIPTOR_FIELD_BASE<0,  0, 1> { };
struct DJDALPHA             : public DESCRIPTOR_FIELD_BASE<0,  1, 0> { };
struct AV_SHEAR             : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct TAU_W                : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct SCALAR               : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct SMAGO_CONST          : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct EFFECTIVE_OMEGA      : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct VELO_GRAD            : public DESCRIPTOR_FIELD_BASE<0,  3, 0> { };
struct FIL_RHO              : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct LOCAL_FIL_VEL_X      : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct LOCAL_FIL_VEL_Y      : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct LOCAL_FIL_VEL_Z      : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct LOCAL_AV_DISS        : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct LOCAL_AV_TKE         : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct LOCAL_SIGMA_ADM      : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct LOCAL_NU_EDDY        : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct FILTERED_VEL_GRAD    : public DESCRIPTOR_FIELD_BASE<0,  3, 0> { };
struct ERROR_COVARIANCE     : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct VARIANCE             : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct TAU_SGS              : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct FILTERED_POPULATION  : public DESCRIPTOR_FIELD_BASE<0,  0, 1> { };
struct INDICATE             : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct BIOGAS_INSTANT       : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct BIOGAS_CUMULATIVE    : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct METHANE_INSTANT      : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct METHANE_CUMULATIVE   : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct CO2_INSTANT          : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };
struct CO2_CUMULATIVE       : public DESCRIPTOR_FIELD_BASE<1,  0, 0> { };

struct TENSOR {
  TENSOR() = delete;

  template <unsigned D, unsigned Q>
  static constexpr unsigned size()
  {
    return (D * (D+1)) / 2; // see `TensorVal` in `core/util.h`
  }

  template <unsigned D, unsigned Q>
  static constexpr unsigned getLocalIndex(const unsigned localIndex)
  {
    return localIndex < size<D,Q>() ? localIndex : throw (std::out_of_range("Index exceeds data field"));
  }
};

//@}

/// \defgroup descriptor_details Descriptor field handling
/// \ingroup descriptor
//@{

template <
  unsigned D,
  unsigned Q,
  typename WANTED_FIELD,
  typename CURRENT_FIELD,
  typename... FIELDS,
  // WANTED_FIELD equals the head of our field list, terminate recursion
  utilities::meta::enable_if_t<std::is_same<WANTED_FIELD,CURRENT_FIELD>::value, int> = 0
>
constexpr unsigned getIndexFromFieldList()
{
  return 0;
}

template <
  unsigned D,
  unsigned Q,
  typename WANTED_FIELD,
  typename CURRENT_FIELD,
  typename... FIELDS,
  // WANTED_FIELD doesn't equal the head of our field list
  utilities::meta::enable_if_t<!std::is_same<WANTED_FIELD,CURRENT_FIELD>::value, int> = 0
>
constexpr unsigned getIndexFromFieldList()
{
  // Break compilation when WANTED_FIELD is not provided by list of fields
  static_assert(sizeof...(FIELDS) > 0, "Field not found.");

  // Add size of current field to implicit offset and continue search
  // for WANTED_FIELD in the tail of our field list
  return CURRENT_FIELD::template size<D,Q>() + getIndexFromFieldList<D,Q,WANTED_FIELD,FIELDS...>();
}


template <unsigned D, unsigned Q>
constexpr unsigned getFieldListSize()
{
  // Field-less descriptor base case
  return 0;
}

template <unsigned D, unsigned Q, typename CURRENT_FIELD, typename... FIELDS>
constexpr unsigned getFieldListSize()
{
  // Calculate size of CURRENT_FIELD and add it to the sum of all remaining field sizes
  return CURRENT_FIELD::template size<D,Q>() + getFieldListSize<D,Q,FIELDS...>();
}

//@}

//@}

}

}

#endif
