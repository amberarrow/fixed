# Fixed -- A simple templated class for fixed-point arithmetic

## Overview

Fixed is a simple templated class for fixed-point arithmetic. On some systems, typically
embedded systems, a floating point unit is either unavailable or is very slow;
fixed-point arithmetic can be much faster on such systems. Additionally, integer adds
and subtracts are often faster than their floating point analogues so if your
computation involves very few multiplies and divides you may find that fixed point
arithmetic yields faster results.

## Usage

A sample program named <b><strong>test_fixed</b></strong> and a *Makefile* are provided.
Just type <b><strong>make</b></strong> to build it and run it with
<b><strong>./test_fixed</b></strong>.

## Contact

I appreciate feedback, so if you have any questions, suggestions or patches, please
send me email: amberarrow at gmail with the usual extension. Thanks.
