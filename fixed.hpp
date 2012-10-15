// A simple templated class for fixed point arithmetic
//
// Fixed point arithmetic is useful in situations where a floating point is
// either not available or is very slow which is often the case for embedded
// systems.
//
// Author: Munagala V. Ramanath
//

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
// for exit()
#include <cstdlib>
#include <math.h>
#include <stdint.h>
#include <type_traits>

using namespace std;

void Fatal( string msg ) __attribute__ ((noreturn));
void
Fatal( string msg )
{
    cerr << msg << endl; exit( 1 );
}

// forward declarations
template< class T, int NFRAC > class Fixed;

template< class T, int NFRAC > ostream &
operator<<( ostream &os, const Fixed< T, NFRAC > &a );

// Multiply without using floats
template< class Fix > Fix
Mul( const Fix &a1, const Fix &a2 )
{
    Fix d = 0;    // return value

    // Almost all the work is done using unsigned arithmetic
    typedef typename Fix::ubase_t Ubase_t;
    register Ubase_t
        uLo, uHi, vLo, vHi, v01, pLo, pHi, a;

    typename Fix::base_t
        u = a1.val, v = a2.val;

    // dispose of easy cases first
    if ( 0 == u || 0 == v) {
        return d;
    }

    // remember if result should be negative and make both operands
    // positive
    //
    bool neg = false;
    if ( u < 0 ) {
        u = -u;
        if ( v > 0 ) {
            neg = true;
        } else {
            v = -v;
        }
    } else if ( v < 0 ) {
        v = -v;
        neg = true;
    }

    // more easy cases
    if ( u == Fix::one ) {
        d.val = neg ? -v : v;
        return d;
    }
    if ( v == Fix::one ) {
        d.val = neg ? -u : u;
        return d;
    }

    // must have v > 0 and u > 0 here
    Ubase_t uz = u, vz = v;
    vLo = vz & Fix::LOMASK,  vHi = vz >> Fix::HALF;
    uLo = uz & Fix::LOMASK,  uHi = uz >> Fix::HALF;

    // If any piece is 0 we save a couple of multiplies
    //
    // Whether these checks are worthwhile depends on the problem
    // domain as well as how expensive integer multiplies and branches
    // are on the particular CPU; with branch prediction, speculative
    // execution, cache effects, etc. it is very difficult to know
    // a priori.
    //
    // If using 32-bit values with an 8-bit fraction, there is a good
    // chance that the high half of one or both operands is 0 if the
    // numbers are small.
    //
    if ( 0 == uHi ) {
        if ( 0 == vHi ) {
            pLo = uLo * vLo; pHi = 0;
        } else if ( 0 == vLo ) {
            a = vHi * uLo; pLo = a << Fix::HALF; pHi = a >> Fix::HALF;
        } else {
            a = vHi * uLo; pHi = a >> Fix::HALF;
            pLo = (vLo * uLo) + (a <<= Fix::HALF);
            if ( pLo < a )
                pHi++;
        }
    } else if ( 0 == vHi ) {
        if ( 0 == uLo ) {
            a = uHi * vLo; pLo = a << Fix::HALF; pHi = a >> Fix::HALF;
        } else {
            a = uHi * vLo; pHi = a >> Fix::HALF;
            pLo = (vLo * uLo) + (a <<= Fix::HALF);
            if ( pLo < a )
                pHi++;
        }
    } else if ( 0 == uLo ) {
        if ( 0 == vLo ) {
            pHi = uHi * vHi; pLo = 0;
        } else {
            pLo = uHi * vLo;
            pHi = (uHi * vHi) + (pLo >> Fix::HALF);
            pLo <<= Fix::HALF;
        }
    } else if ( 0 == vLo ) {
        pLo = vHi * uLo;
        pHi = (uHi * vHi) + (pLo >> Fix::HALF);
        pLo <<= Fix::HALF;
    } else {

        /* All pieces are non-zero */
        pLo = uLo * vLo;  /* low half of product */
        pHi = uHi * vHi;  /* high half of product */

        // compute middle product p01
        if ( vHi < vLo ) {

            v01 = vLo - vHi;
            if ( uHi > uLo ) {
                /* p01 > 0 */
                a = pHi + (uHi - uLo) * v01;
                if ( a < pHi )
                    pHi += Fix::CARRY;
                a += pLo;
                if ( a < pLo )
                    pHi += Fix::CARRY;
            } else {
                /* p01 <= 0, |p01| < p0 */
                a = pLo - (uLo - uHi) * v01;
                a += pHi;
                if ( a < pHi )
                    pHi += Fix::CARRY;
            }

        } else {                /* vHi >= vLo */

            v01 = vHi - vLo;
            if ( uHi > uLo ) {
                /* p01 <= 0, |p01| < pHi */
                a = pHi - (uHi - uLo) * v01;
                a += pLo;
                if ( a < pLo )
                    pHi += Fix::CARRY;
            } else {
                /* p01 >= 0 */
                a = pHi + (uLo - uHi) * v01;
                if ( a < pHi )
                    pHi += Fix::CARRY;
                a += pLo;
                if ( a < pLo )
                    pHi += Fix::CARRY;
            }
        }

        // add in the middle product
        pHi += (a >> Fix::HALF);
        pLo += (a <<= Fix::HALF);
        if ( pLo < a )
            pHi++;
    }

    // round result to NFRAC fraction bits
    pLo += Fix::half;
    if ( pLo < static_cast<Ubase_t>(Fix::half) )
        pHi++;

    // if the high bits that will be discarded are non-zero, we have
    // overflow; likewise if the high-bit of the result is non-zero
    //
    if ( pHi & (Fix::intMask | Fix::half) ) {    // overflow
        cerr << "pHi = " << pHi << " [0x" << hex << pHi << "], pLo = "
             << pLo << dec << endl;
        cerr << "u = "  << u << ", v = " << v << ", d.val = " << d.val
             << endl;
        Fatal( "ERROR(Fixed::Mul): overflow" );
    }

    // deal with the sign
    d.val = (pLo >> Fix::N)
            | ((pHi & Fix::fracMask) << (Fix::NBITS - Fix::N));
    if ( neg ) {
        d.val = -d.val;
    }

    return d;
}                                                                 // Mul

//---------------------------------------------------------------  Fixed
//
// A simple fixed-point arithmetic class.
//
// First argument T is the signed integer type to use (usually 'int',
// 'long', or 'long long'
//
// Second argument NFRAC is the number of fraction bits.
//
// Must have 1 <= NFRAC <= (nBits(T) - 2) where nBits(T) is the number
// of bits in T. For example, if T is 'int', and an int has 32 bits, we
// must have 1 <= NFRAC <= 30
//
template< class T, int NFRAC >
class Fixed {

    friend ostream & operator<< < T, NFRAC > (
        ostream &os, const Fixed< T, NFRAC > &a );

    friend Fixed Mul< Fixed > (
        const Fixed &a1, const Fixed &a2 );

public:

    typedef T base_t;
    typedef typename std::make_unsigned<T>::type ubase_t;

    // Some useful constants, smallest and largest values,
    // mask for fraction, etc.
    //
    static const base_t
        one, half, negOne, negHalf, min, max, fracMask, intMask;
    static const ubase_t
        CARRY, LOMASK;
    static const int
        N     = NFRAC,             // no. of fraction bits
        NBITS = sizeof( T ) * 8,   // no. of bits
        HALF  = NBITS >> 1 ;       // half the no. of bits

    // fmult  --  the multiplier for converting a double to a Fixed; if
    //            all bits are fraction bits, (1 << NFRAC) will be zero,
    //            so we shift by 1 less and multiply by 2 as a double
    // fmaxS  --  max as a double (true maximum scaled by fmult)
    // fminS  --  min as a double (true minimum scaled by fmult)
    // fmin   --  min as a double scaled down by fmult (true minimum)
    // fmax   --  max as a double scaled down by fmult (true maximum)
    //
    static const double
        fmult, fmaxS, fminS, fmax, fmin;

    // ctors
    Fixed( ) { }
    Fixed( const Fixed &a ) : val( a.val ) { }

    Fixed( const double &a ) {
#ifdef DEBUG
        if ( a > fmax || a < fmin ) {
            Fatal( "ERROR(Fixed(double)): Overflow" );
        }
#endif
        val = static_cast<T>( a * fmult );
    }

    Fixed( const int &a ) : val( a << NFRAC ) {
#ifdef DEBUG
        if ( (val >> NFRAC) != a ) {
            Fatal( "ERROR(Fixed(int)): Overflow" );
        }
#endif
    }
    Fixed( const long &a ) : val( a << NFRAC ) {
#ifdef DEBUG
        if ( (val >> NFRAC) != a ) {
            Fatal( "ERROR(Fixed(long)): Overflow" );
        }
#endif
    }

    // No dtor needed

    // operators
    bool operator<(  const Fixed &a ) { return val < a.val; }
    bool operator<=( const Fixed &a ) { return val <= a.val; }
    bool operator>(  const Fixed &a ) { return val > a.val; }
    bool operator>=( const Fixed &a ) { return val >= a.val; }
    bool operator==( const Fixed &a ) { return val == a.val; }
    bool operator!=( const Fixed &a ) { return val != a.val; }

    Fixed operator<<( int a ) { return val << a; }
    Fixed operator>>( int a ) { return val >> a; }

    Fixed operator=( const Fixed &a ) { val = a.val; return *this; }

    Fixed operator+=( const Fixed &a ) {
#ifdef DEBUG
        T tmp = val + a.val;
        if ( (tmp < 0 && val > 0 && a.val > 0)
                || (tmp > 0 && val < 0 && a.val < 0) ) {
            Fatal( "ERROR(Fixed::[+= or +]): Overflow" );
        }
#endif
        val += a.val;
        return *this;
    }

    Fixed operator-=( const Fixed &a ) {
#ifdef DEBUG
        T tmp = val - a.val;
        if ( (tmp < 0 && val > 0 && a.val < 0)
                || (tmp > 0 && val < 0 && a.val > 0) ) {
            Fatal( "ERROR(Fixed::[-= or -]): Overflow" );
        }
#endif
        val -= a.val;
        return *this;
    }

    // conversion to double is a quick and dirty implementation which
    // may lose some fraction bits if T is a 64-bit value; more exact
    // all-integer methods exist (do later)
    //
    Fixed operator*=( const Fixed &a ) {
        double tmp = (val / fmult) * a.val;
#ifdef DEBUG
        if ( tmp > fmaxS || tmp < fminS ) {
            Fatal( "ERROR(Fixed::[*= or *]): Overflow" );
        }
#endif
        val = static_cast< T >( tmp );
        return *this;
    }

    // conversion to double is a quick and dirty implementation which
    // may lose some fraction bits if T is a 64-bit value;; more exact
    // all-integer methods exist (do later)
    //
    Fixed operator/=( const Fixed &a ) {
        double tmp = static_cast< double >( val ) / a.val;
#ifdef DEBUG
        if ( tmp > fmax || tmp < fmin ) {
            Fatal( "ERROR(Fixed::[/= or /]): Overflow" );
        }
#endif
        val = static_cast< T >( tmp * fmult );
        return *this;
    }

    operator double() { return Fval(); };

    // other methods

    // test if integer
    bool IsInt() const { return val & fracMask ? false : true; }

    // integer part
    T Int() const { return val >> NFRAC; }

    // NOTE: many of these are destructive

    // Round to nearest integer.
    //
    // Note that simply zero-ing out the fraction bits corresponds to
    // truncation towards -infinity, which makes rounding easy here
    //
    void Round() {
        if ( 0 == (val & fracMask) ) {    // no fraction
            return;
        }
        val += half;
        val &= intMask;
    }

    // RoundDown (nearest integer toward 0).
    //
    // Note that simply zero-ing out the fraction bits corresponds to
    // truncation towards -infinity, so negative values need special
    // treatment.
    //
    void RoundDown() {
        if ( 0 == (val & fracMask) ) {    // no fraction
            return;
        }
        if ( val < 0 ) {
            val += one;
        }
        val &= intMask;
    }

    // RoundUp (nearest integer away from 0).
    //
    // Note that simply zero-ing out the fraction bits corresponds to truncation
    // towards -infinity, so positive values need special treatment.
    //
    void RoundUp() {
        if ( 0 == (val & fracMask) ) {    // no fraction
            return;
        }
        if ( val > 0 ) {
            val += one;
        }
        val &= intMask;
    }

    // Floor (nearest integer toward -infinity).
    //
    void Floor() {
        val &= intMask;
    }

    // Ceiling (nearest integer toward +infinity).
    //
    void Ceil() {
        if ( 0 == (val & fracMask) ) {    // no fraction
            return;
        }
        val += one;
        val &= intMask;
    }

    // value as a double (some loss of precision possible if T is 64-bit)
    double Fval() const {
        return val / fmult;
    }

    // square root
    Fixed Sqrt() {
        return fmult * sqrt( static_cast< double >( val ) / fmult );
    }

private:
    T val;

};  // Fixed

// Definitions of constants
//
// Note that NFRAC is required to be at most (nbits - 2) where nbits is
// the number of bits in the base type, so (1 << NFRAC) will never be
// zero or negative.
//
template< class T, int NFRAC > const T
    Fixed<T,NFRAC>::one      = static_cast< T >( 1 ) << NFRAC;
template< class T, int NFRAC > const T
    Fixed<T,NFRAC>::half     = Fixed<T,NFRAC>::one >> 1;
template< class T, int NFRAC > const T
    Fixed<T,NFRAC>::negOne   = static_cast< T >( -1 ) << NFRAC;
template< class T, int NFRAC > const T
    Fixed<T,NFRAC>::negHalf  = Fixed<T,NFRAC>::negOne >> 1;
template< class T, int NFRAC > const T
    Fixed<T,NFRAC>::min      = static_cast< T >( 1 )
                                  << ((sizeof( T ) << 3) - 1);
template< class T, int NFRAC > const T
    Fixed<T,NFRAC>::max      = Fixed<T,NFRAC>::min ^ static_cast< T >( -1 );

template< class T, int NFRAC > const T
    Fixed<T,NFRAC>::fracMask = Fixed<T,NFRAC>::one - 1;
template< class T, int NFRAC > const T
    Fixed<T,NFRAC>::intMask  = ~Fixed<T,NFRAC>::fracMask;

template< class T, int NFRAC > const typename Fixed<T,NFRAC>::ubase_t
    Fixed<T,NFRAC>::CARRY = static_cast< Fixed<T,NFRAC>::ubase_t >( 1 ) << HALF;
template< class T, int NFRAC > const typename Fixed<T,NFRAC>::ubase_t
    Fixed<T,NFRAC>::LOMASK = Fixed<T,NFRAC>::CARRY - 1;

template< class T, int NFRAC > const double
    Fixed<T,NFRAC>::fmult  = Fixed<T,NFRAC>::one;
template< class T, int NFRAC > const double
    Fixed<T,NFRAC>::fmaxS  = Fixed<T,NFRAC>::max;
template< class T, int NFRAC > const double
    Fixed<T,NFRAC>::fminS  = Fixed<T,NFRAC>::min;
template< class T, int NFRAC > const double
    Fixed<T,NFRAC>::fmax   = Fixed<T,NFRAC>::max
                                / Fixed<T,NFRAC>::fmult;
template< class T, int NFRAC > const double
    Fixed<T,NFRAC>::fmin   = Fixed<T,NFRAC>::min
                                / Fixed<T,NFRAC>::fmult;

// function definitions

template< class T, int NFRAC > inline Fixed< T, NFRAC >
operator+( const Fixed< T, NFRAC > &a1,
           const Fixed< T, NFRAC > &a2 )
{
    Fixed< T, NFRAC > r( a1 );
    r += a2;
    return r;
}  // operator+

template< class T, int NFRAC > inline Fixed< T, NFRAC >
operator-( const Fixed< T, NFRAC > &a1,
           const Fixed< T, NFRAC > &a2 )
{
    Fixed< T, NFRAC > r = a1;
    r -= a2;
    return r;
}  // operator-

template< class T, int NFRAC > inline Fixed< T, NFRAC >
operator*( const Fixed< T, NFRAC > &a1,
           const Fixed< T, NFRAC > &a2 )
{
    Fixed< T, NFRAC > r = a1;
    r *= a2;
    return r;
}  // operator*

template< class T, int NFRAC > inline Fixed< T, NFRAC >
operator/( const Fixed< T, NFRAC > &a1,
           const Fixed< T, NFRAC > &a2 )
{
    Fixed< T, NFRAC > r = a1;
    r /= a2;
    return r;
}  // operator/

// output a fixed point number (some loss of precision possible)
template< class T, int NFRAC > ostream &
operator<<( ostream &os, const Fixed< T, NFRAC > &a )
{
    return os << a.Fval();
}  // operator<<

// instantiations
typedef Fixed< int64_t, 8 > Fix8;
typedef Fixed< int64_t, 16 > Fix16;
