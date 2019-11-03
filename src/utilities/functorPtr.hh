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

#ifndef FUNCTOR_PTR_HH
#define FUNCTOR_PTR_HH

#include "functorPtr.h"

#include "io/ostreamManager.h"

namespace olb {

template<typename F>
FunctorPtr<F>::FunctorPtr(F& f):
  _ownF(),
  _sharedF(),
  _f(&f),
  _owning(false)
{ }

template<typename F>
FunctorPtr<F>::FunctorPtr(F* f):
  _ownF(),
  _sharedF(),
  _f(f),
  _owning(false)
{ }

template<typename F>
FunctorPtr<F>::FunctorPtr(const std::unique_ptr<F>& f):
  FunctorPtr(f.get())
{ }

template<typename F>
FunctorPtr<F>::FunctorPtr(std::unique_ptr<F>&& f):
  _ownF(std::move(f)),
  _f(_ownF.get()),
  _owning(true)
{ }

template<typename F>
FunctorPtr<F>::FunctorPtr(std::shared_ptr<F> f):
  _sharedF(std::move(f)),
  _f(_sharedF.get()),
  _owning(true)
{ }

template<typename F>
FunctorPtr<F>::FunctorPtr():
  FunctorPtr(nullptr)
{ }

template<typename F>
typename std::add_lvalue_reference<F>::type FunctorPtr<F>::operator*() const
{
  return *_f;
}

template<typename F>
typename std::add_pointer<F>::type FunctorPtr<F>::operator->() const noexcept
{
  return _f;
}

template<typename F>
FunctorPtr<F>::operator bool() const
{
  return _f != nullptr;
}

template<typename F>
bool FunctorPtr<F>::isOwning() const
{
  return _owning;
}

template<typename F>
template<typename G>
auto FunctorPtr<F>::toShared() -> typename std::enable_if< util::has_identity_functor<G>::value, std::shared_ptr<F>>::type
{
  if ( _owning && _ownF ) {
    // convert unique_ptr to shared_ptr
    _sharedF = std::shared_ptr<F>(std::move(_ownF));
    _f       = _sharedF.get();
    return _sharedF;
  }
  else if ( _sharedF ) {
    // either we are owning or toShared was called previously and
    // prepared a F::identity_functor_type we can reuse here:
    return _sharedF;
  }
  else {
#ifdef OLB_DEBUG_UNSAFE_FUNCTOR_ARITHMETIC
    OstreamManager clerr(std::cerr,"FunctorPtr");
    clerr << "WARNING: Non-owning functor wrapped in std::shared_ptr. Lifetime of composed functor is not guaranteed." << std::endl;
#endif

    // wrap non-owning pointer in identity functor
    _sharedF = std::shared_ptr<F>(new typename F::identity_functor_type(*_f));
    return _sharedF;
  }
}

template<typename F>
template<typename G>
auto FunctorPtr<F>::toShared() -> typename std::enable_if<!util::has_identity_functor<G>::value, std::shared_ptr<F>>::type
{
  if ( _owning && _ownF ) {
    // convert unique_ptr to shared_ptr
    _sharedF = std::shared_ptr<F>(std::move(_ownF));
    _f       = _sharedF.get();
    return _sharedF;
  }
  else if ( _sharedF ) {
    return _sharedF;
  }
  else {
    OstreamManager clerr(std::cerr,"FunctorPtr");
    clerr << "ERROR: No identity functor defined. Required to lift functor reference into std::shared_ptr arithmetic." << std::endl;
    std::exit(-1);
  }
}

}

#endif
