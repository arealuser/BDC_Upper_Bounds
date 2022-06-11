SOURCES = $(wildcard *.cc)
HEADERS = $(wildcard *.h)

MAIN_SOURCES = $(wildcard *.cpp)
MAINS = $(MAIN_SOURCES:%.cpp=%.out)

OBJECTS = $(SOURCES:%.cc=%.o)

# CC := $(shell which clang || which gcc)
CC := g++-9
CFLAGS = -Wall -W -O3 -fno-exceptions -fno-rtti -std=c++2a
LIBS = stdc++ m
LDFLAGS = $(LIBS:%=-l%)

bit_channel: $(OBJECTS) all_mains
	echo done

all_mains : $(OBJECTS) $(MAINS)
		echo done

test: all_mains
	./test_bit_baa.out
	./test_baa.out
	./test_transition_probability_computation.out

%.out : %.cpp $(OBJECTS)
		$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $^

%.o : %.cc $(HEADERS)
		$(CC) $(CFLAGS) -c -o $@ $<

.PHONY : clean
clean :
		rm -f $(OBJECTS) $(MAINS)