/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef FUNCTOR_PTR_H
#define FUNCTOR_PTR_H

#include <memory>
#include <type_traits>

namespace olb {

namespace util {

/// Backport of C++17's std::void_t
template<typename...>
struct void_t {
  using type = void;
};

/// Indicates existence of F::identity_functor_type typedef
template<class F, class U = void>
struct has_identity_functor : std::false_type { };
/// Indicates existence of F::identity_functor_type typedef
template<class F>
struct has_identity_functor<F, typename void_t<typename F::identity_functor_type>::type> : std::true_type { };

}

/// Smart pointer for managing the various ways of passing functors around
/**
 * There is a rich set of functors that compose other functors by accepting
 * them as constructor arguments. e.g. SuperLpNorm3D, SuperPlaneIntegralF3D
 *
 * Previously all of these functors only accepted other functors by reference
 * which prevented e.g. the use case of accepting a std::unique_ptr to a
 * indicator freshly created by SuperGeometry3D::getMaterialIndicator.
 *
 * This class hides the memory management required to support accepting either
 * functor references, non-owning pointers or owning pointers in a single argument.
 *
 * Note that if this class is constructed by reference or raw pointer the caller
 * is responsible for assuring that the passed functor is available for the full
 * lifetime. This is not the case when constructed by a smart pointer.
 **/
template<typename F>
class FunctorPtr {
private:
  /// Optional functor owner
  std::unique_ptr<F> _ownF;
  /// Optional shared functor
  std::shared_ptr<F> _sharedF;
  /// Pointer to the exposed functor
  F* _f;
  /// True iff FunctorPtr was constructed owning
  const bool _owning;

public:
  /// Constructor for transparently accepting a functor reference
  FunctorPtr(F& f);
  /// Constructor for transparently accepting a non-owning functor pointer
  FunctorPtr(F* f);
  /// Constructor for transparently accepting a owning functor pointer
  /**
   * Equivalent to construction using a non-owning functor pointer
   **/
  FunctorPtr(const std::unique_ptr<F>& f);
  /// Constructor for transparently taking ownership of a owning functor pointer
  FunctorPtr(std::unique_ptr<F>&& f);
  /// Constructor for transparently sharing ownership of a functor
  FunctorPtr(std::shared_ptr<F> f);
  /// Constructor for an empty instance
  /**
   * Equivalent to construction by nullptr.
   **/
  FunctorPtr();

  /// Copy construction is disabled as it is incompatible with unique ownership
  FunctorPtr(FunctorPtr&) = delete;
  /// Move constructor
  FunctorPtr(FunctorPtr&&) = default;

  /// Perfect forwarding functor operator
  template<typename... Args>
  bool operator()(Args... args)
  {
    return _f->operator()(std::forward<Args>(args)...);
  }

  /// \returns reference to the exposed functor
  typename std::add_lvalue_reference<F>::type operator*() const;
  /// Enable pointer-like access to the exposed functor's members
  typename std::add_pointer<F>::type operator->() const noexcept;

  /// Indicates whether a functor instance is exposed
  /**
   * \returns true iff a functor is exposed
   *
   * This is useful for supporting optional functors in constructor interfaces
   * by constructing FunctorPtr using a nullptr
   **/
  operator bool() const;

  /// Indicates whether a functor instance is both exposed and owned
  /**
   * \returns true iff FunctorPtr took ownership of managed functor during construction
   *
   * It is not enough to check _ownF and _sharedF for non-nullptr contents as _sharedF
   * is also used to store F::identity_functor_type instances created by toShared.
   **/
  bool isOwning() const;

  /// Lifts managed functor reference for usage in std::shared_ptr<F> arithmetic
  template<typename G = F>
  auto toShared() -> typename std::enable_if< util::has_identity_functor<G>::value, std::shared_ptr<F>>::type;
  /// Lifts managed functor reference for usage in std::shared_ptr<F> arithmetic (limited)
  template<typename G = F>
  auto toShared() -> typename std::enable_if<!util::has_identity_functor<G>::value, std::shared_ptr<F>>::type;

};

}

#endif
