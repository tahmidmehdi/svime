/*
 * distribution.h
 * Statistical Distributions
 *  Updated on: Feb 21, 2019
 *      Author: Tahmid Mehdi
 *
Copyright 2019 Tahmid Mehdi
This file is part of svime.

svime is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

svime is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with svime.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_DISTRIBUTION_H_
#define SRC_DISTRIBUTION_H_
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Dense>

namespace Dist {
// PDF of Beta(x|a, b)
float beta_pdf(const float x, const float a, const float b) {
	return std::pow(x, a-1)*std::pow(1-x, b-1)/boost::math::beta(a, b);
}
// log PDF of Beta(x|a, b). lnx is log(x) and lnix is log(1-x)
float ln_beta(const double lnx, const double lnix, const double a, const double b) {
	return (a-1)*lnx+(b-1)*lnix-(boost::math::lgamma(a)+boost::math::lgamma(b)-boost::math::lgamma(a+b));
}
// PDF of Multinomial(x|p) where x is a base and p contains base probabilities
float multinomial_pmf(const char x, Eigen::Vector4d p) {
	return std::pow(p(0), (x=='A'))*std::pow(p(1), (x=='C'))*std::pow(p(2), (x=='G'))*std::pow(p(3), (x=='T'));
}
// log PDF of Multinomial(x|p). lnp contains log base probabilities
float ln_multinomial(const char x, Eigen::Vector4d lnp) {
	return (x=='A')*lnp(0)+(x=='C')*lnp(1)+(x=='G')*lnp(2)+(x=='T')*lnp(3);
}
// PDF of Dirichlet(x|a) where a is vector of concentrations
float dirichlet_pdf(Eigen::Vector4d x, Eigen::Vector4d a) {
	Eigen::Vector4d aGamma;
	for ( int i = 0; i < 4; ++i ) {
		aGamma(i) = std::pow(x(i), a(i)-1)/boost::math::tgamma(a(i));
	}
	return boost::math::tgamma(a.sum())*aGamma.prod();
}
// log PDF of Dirichlet(p|a). lnp contains log base probabilities
float ln_dirichlet(Eigen::Vector4d lnp, Eigen::Vector4d a) {
	Eigen::Vector4d aGamma;
	for ( int i = 0; i < 4; ++i ) {
		aGamma(i) = (a(i)-1)*lnp(i)-boost::math::lgamma(a(i));
	}
	return boost::math::lgamma(a.sum())+aGamma.sum();
}
}

#endif /* SRC_DISTRIBUTION_H_ */
