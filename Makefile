CXXFLAGS = -g -O2 -Wall -Wextra -std=c++0x

test_fixed: test_fixed.cpp

clean:
	rm -f test_fixed *.o
