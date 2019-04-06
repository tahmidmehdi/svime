//****************************************************************************80
//
//  Purpose:
//
//    DIGAMMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 March 2016
//
//  Author:
//
//    Original FORTRAN77 version by Jose Bernardo.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jose Bernardo,
//    Algorithm AS 103:
//    Psi ( Digamma ) Function,
//    Applied Statistics,
//    Volume 25, Number 3, 1976, pages 315-317.
//
#ifndef SRC_ASA103_HPP_
#define SRC_ASA103_HPP_

double digamma ( double x, int *ifault );
void psi_values ( int *n_data, double *x, double *fx );
void timestamp ( );

#endif /* SRC_ASA103_HPP_ */
