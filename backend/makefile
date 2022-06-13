SOURCES = $(wildcard *.cc)
HEADERS = $(wildcard *.h)

TRANSITION_SOURCES = $(wildcard transition_counts/*.cc)
TRANSITION_HEADERS = $(wildcard transition_counts/*.h) $(wildcard transition_counts/*.inl)
TRANSITION_OBJ = transition_cache.o

MAIN_SOURCES = $(wildcard *.cpp)
MAINS = $(MAIN_SOURCES:%.cpp=%.out)

OBJECTS = $(SOURCES:%.cc=%.o)

# CC := $(shell which clang || which gcc)
CC := g++-9
CFLAGS = -Wall -W -O3 -fno-exceptions -fno-rtti -std=c++2a -g
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

%.out : %.cpp $(OBJECTS) $(TRANSITION_OBJ)
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $(OBJECTS) $(TRANSITION_OBJ) $<
%.o : %.cc $(HEADERS) 		
	$(CC) $(CFLAGS) -c -o $@ $<

$(TRANSITION_OBJ): $(TRANSITION_SOURCES) $(TRANSITION_HEADERS)
		$(CC) $(CFLAGS) -c -o $(TRANSITION_OBJ) transition_counts/cached_transition_count.cc

.PHONY : clean
clean :
	rm -f $(OBJECTS) $(MAINS) $(TRANSITION_OBJ)