/*
 * transport_equations.h
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#pragma once

#include <boost/version.hpp>
#if BOOST_VERSION >= 107200
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#else
#include <boost/math/interpolators/cubic_b_spline.hpp>
#endif

namespace BSMPT {
namespace Baryo {
#if BOOST_VERSION >= 107200
template<typename T>
using boost_cubic_b_spline = boost::math::interpolators::cardinal_cubic_b_spline<T>;
#else
template<typename T>
using boost_cubic_b_spline = boost::math::cubic_b_spline<T>;
#endif
}
}
