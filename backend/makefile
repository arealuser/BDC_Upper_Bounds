SOURCES = $(wildcard *.cc)
HEADERS = $(wildcard *.h)

MAIN_SOURCES = $(wildcard *.cpp)
MAINS = $(MAIN_SOURCES:%.cpp=%.out)

OBJECTS = $(SOURCES:%.cc=%.o)

# CC := $(shell which clang || which gcc)
CC := g++-9
CFLAGS = -Wall -W -O3 -fno-exceptions -fno-rtti -std=c++2a -g -D__SOURCE_PATH__="\"$$(pwd)\"" -fsanitize=address
LIBS = stdc++ m
LDFLAGS = $(LIBS:%=-l%)

bit_channel: all_mains
	echo done

all_mains : $(MAINS)
	echo done

all: $(MAINS) $(OBJECTS)
	echo done

test: test_bit_baa.out test_transition_probability_computation.out test_baa.out
	./test_transition_probability_computation.out
	./test_bit_baa.out
	./test_baa.out	

%.out : %.cpp $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $(OBJECTS) $<

%.o : %.cc $(HEADERS) 		
	$(CC) $(CFLAGS) -c -o $@ $<



.PHONY : clean
clean :
	rm -f $(OBJECTS) $(MAINS)