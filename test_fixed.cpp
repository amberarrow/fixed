// Tests for fixed point arithmetic
//
// Author: Munagala V. Ramanath
//

#include <sys/time.h>
#include <ctime>

#include "fixed.hpp"

// instantiations
typedef Fixed< int64_t, 8 > Fix8;
typedef Fixed< int64_t, 16 > Fix16;

static struct timeval tv;    // scratch

struct Time {    // simple time class for time trials
    uint64_t val;

    Time() {    // initialize with current time
        gettimeofday( &tv, NULL );
        val = tv.tv_usec  + tv.tv_sec * 1000000;
    }

    // print difference between current and saved time in milliseconds
    void diff( char const *const msg ) const {
        Time cur;
        cout << msg << (cur.val - val) / 1000 << " ms" << endl;
    }  // diff
};  // Time

// return random float in [-1000,1000]
inline double frand() {
    return -1000.0 + 2000.0 * random() / (double)RAND_MAX;
}  // frand

// test
int
main( int argc, char **argv )
{
//#define FIX Fix8
#define FIX Fix16

    FIX f1( 10 ), f2( 20 ), f3 = f1 * f2, f4 = f1 / f2, f5 = f1 + f2,
        f6 = f1 - f2;

    cout << "max = " << FIX::max << ", min = " << FIX::min << endl
         << "fmax = " << FIX::fmax
         << ", fmin = " << FIX::fmin << endl
         << "sizeof( long ) = " << sizeof(long) << endl
         << "fmult = " << FIX::fmult
         << ", f1 = " << f1 << ", f2 = " << f2
         << ", f3 = " << f3 << ", f4 = " << f4
         << ", f5 = " << f5 << ", f6 = " << f6
         << endl;

    // time trial
    Time t1;
    int const n = argc > 1 ? atoi( argv[ 1 ] ) : 1000000;
    double sum = 0.0;    // prevent gcc optimizing away whole loop
    for ( int i = 0; i < n; ++i ) {    // use floats to multiply
        f1 = frand(), f2 = frand(), f3 = f1 * f2;
        sum += f3;
    }
    t1.diff( "With floats: " );
    for ( int i = 0; i < n; ++i ) {    // no floats
        f1 = frand(), f2 = frand(), f3 = Mul( f1, f2 );
        sum += f3;
    }
    t1.diff( "No floats: " );
    cout << "sum = " << sum << endl;

    return 0;
}  // main
