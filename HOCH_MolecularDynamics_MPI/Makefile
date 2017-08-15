EXECUTABLE := HOCH_MolecularDynamics_MPI

CODEROOT := ..
SRCDIR:= ./src
OBJDIR:= ./build

SHELL = /bin/sh

WARNINGS := -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Werror -Wuninitialized -Wconversion
CFLAGS := -g -O2 -DHAVE_INLINE
CFLAGS += -I$(HOME)/local/include -I$(CODEROOT)/shared/include
LDFLAGS := -lLinkedList -lFormaldehyde -lCalculus 

# Use for GCC
CC := gcc
CFLAGS += $(WARNINGS)
CFLAGS += -std=c99 -fshort-enums -fno-common
LDFLAGS += -static -lgsl -lcblas -lblas -lgfortran -lm -L$(HOME)/local/lib -L$(CODEROOT)/shared/lib
LDFLAGS += -L/soft/libraries/alcf/current/gcc/CBLAS/lib/
LDFLAGS += -L/soft/libraries/alcf/current/gcc/BLAS/lib/
LDFLAGS += -L/soft/apps/ibmcmp-feb2013/xlf/bg/14.1/bglib64 -lxlopt -lxl -lxlf90_r -lxlfmath

# Use for XL
# CC := mpixlc_r
# CFLAGS += -qlanglvl=stdc99
# LDFLAGS += -lgsl -lcblas -lblas -lm -L$(HOME)/cross/lib -L$(CODEROOT)/shared/lib
# LDFLAGS += -L/soft/libraries/alcf/current/xl/CBLAS/lib
# LDFLAGS += -L/soft/libraries/alcf/current/xl/BLAS/lib
# LDFLAGS += -L/soft/apps/ibmcmp-feb2013/xlf/bg/14.1/bglib64 -lxlopt -lxl -lxlf90_r -lxlfmath


HEADERS := $(wildcard $(SRCDIR)/*.h)
SOURCES := $(wildcard $(SRCDIR)/*.c)
OBJECTS := $(patsubst $(SRCDIR)%.c,$(OBJDIR)%.o, $(SOURCES))

all: $(SOURCES) $(HEADERS) $(EXECUTABLE)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) -O2 $< -o $@

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

clean:
	rm $(EXECUTABLE) || true
	cd $(OBJDIR) && rm * || true

install: $(EXECUTABLE)
	install -m 0755 $(EXECUTABLE) $(CODEROOT)/shared/bin

remove:
	rm  $(CODEROOT)/shared/bin/$(EXECUTABLE)
