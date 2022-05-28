SOURCES = $(wildcard *.cc)
HEADERS = $(wildcard *.h)

MAIN_SOURCES = $(wildcard *.cpp)
MAINS = $(MAIN_SOURCES:%.cpp=%.out)

OBJECTS = $(SOURCES:%.cc=%.o)

# CC := $(shell which clang || which gcc)
CC := g++
CFLAGS = -Wall -W -O3 -fno-exceptions -fno-rtti -std=c++2a
LIBS = stdc++ m
LDFLAGS = $(LIBS:%=-l%)

all : $(OBJECTS) $(MAINS)
		echo done


%.out : %.cpp $(OBJECTS)
		$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $^

%.o : %.cc $(HEADERS)
		$(CC) $(CFLAGS) -c -o $@ $<

.PHONY : clean
clean :
		rm -f $(OBJECTS) $(MAINS)