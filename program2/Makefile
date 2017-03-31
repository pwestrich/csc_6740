
#utility information
CXX      := g++
CPPFLAGS := -std=c++11 -O0 -pedantic -Wall -Wextra -g
LDFLAGS  := $(CPPFLAGS)
AR 	 := ar -r
RM 	 := rm -f
MV 	 := mv

#build information
INC_DIRS := -I./include
LIB_DIRS :=
LIBS 	 :=

#how we build things
COMPILE := $(CXX) $(CPPFLAGS) $(INC_DIRS) -c
LINK    := $(CXX) $(LDFLAGS) $(LIB_DIRS) -o

#what to build
TARGET := main

#files
HEADER := $(wildcard *.h)
SOURCE := $(wildcard *.cpp)
OBJECT := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

.PHONY: all clean

all: $(TARGET)

clean:

	$(RM) $(TARGET) $(OBJECT)
	$(RM) -r *.dYSM

main: main.o

	$(LINK) main main.o $(LIBS)

main.o: main.cpp

	$(COMPILE) main.cpp

#generic rule for anything missed
%.o: %.cpp %.h

	$(COMPILE) $<
