CPPFLAGS=-I../
CFLAGS= -g -Wall -O2 -finline-functions
CXXFLAGS= $(CFLAGS)

OBJS=quartic.o

libmathbits.a: $(OBJS)
	$(AR) rvs $@ $^