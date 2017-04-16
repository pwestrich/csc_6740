
#utility information
CXX      := mpicc
CPPFLAGS := -std=c11 -O0 -pedantic -Wall -Wextra -g
LDFLAGS  := $(CPPFLAGS)
AR 	 := ar -r
RM 	 := rm -f
MV 	 := mv

#build information
INC_DIRS :=
LIB_DIRS :=
LIBS 	 :=

#how we build things
COMPILE := $(CXX) $(CPPFLAGS) $(INC_DIRS) -c
LINK    := $(CXX) $(LDFLAGS) $(LIB_DIRS) -o

#what to build
TARGET := main

#files
HEADER := $(wildcard *.h)
SOURCE := $(wildcard *.c)
OBJECT := $(patsubst %.c,%.o,$(wildcard *.c))

.PHONY: all clean

all: $(TARGET)

clean:

	$(RM) $(TARGET) $(OBJECT)
	$(RM) -r *.dYSM

$(TARGET): $(OBJECT)

	$(LINK) $(TARGET) $(OBJECT) $(LIBS)

main.o: main.c getRealTime.h

	$(COMPILE) main.c

#generic rule for anything missed
%.o: %.c %.h

	$(COMPILE) $<
