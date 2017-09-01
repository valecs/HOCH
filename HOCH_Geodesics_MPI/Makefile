EXECUTABLE := HOCH_Geodesics_MPI

CODEROOT := ..
SRCDIR:= ./src
TSTDIR:= ./test
OBJDIR:= ./build

SHELL = /bin/sh

CFLAGS := -g -O2 -DHAVE_INLINE
CFLAGS += -I$(HOME)/toor/include -I$(HOME)/local/include -I$(CODEROOT)/shared/include

WARNINGS := -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Werror -Wuninitialized -Wconversion

LDFLAGS := -lLinkedList -lFormaldehyde -lCalculus  -lgsl -lm -lgfortran
LDFLAGS += -L$(CODEROOT)/shared/lib -L$(HOME)/toor/lib

XLKEY := /wrappers/xl
ifeq ($(XLKEY),$(findstring $(XLKEY), $(PATH)))
	CC := mpixlc_r
	CFLAGS += -qlanglvl=stdc99
	LDFLAGS += -lcblas -lblas
	LDFLAGS += -L$(HOME)/cross/lib
	LDFLAGS += -L/soft/libraries/alcf/current/xl/CBLAS/lib
	LDFLAGS += -L/soft/libraries/alcf/current/xl/BLAS/lib
	LDFLAGS += -L/soft/apps/ibmcmp-feb2013/xlf/bg/14.1/bglib64 -lxlopt -lxl -lxlf90_r -lxlfmath
	LDFLAGS += -s
else # Assume GCC or compatible
	CC := mpicc
	CFLAGS += -std=c99 -fshort-enums -fno-common
	CFLAGS += $(WARNINGS)
	LDFLAGS += -lgslcblas
endif

HEADERS := $(wildcard $(SRCDIR)/*.h)
SOURCES := $(wildcard $(SRCDIR)/*.c)
OBJECTS := $(patsubst $(SRCDIR)%.c,$(OBJDIR)%.o, $(SOURCES))
TSOURCE := $(wildcard $(TSTDIR)/*.c)
TOBJECT := $(patsubst $(TSTDIR)%.c,$(TSTDIR)%.o, $(TSOURCE))

all: $(SOURCES) $(HEADERS) $(EXECUTABLE)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(TSTDIR)/%.o : $(TSTDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

clean:
	rm $(EXECUTABLE) || true
	cd $(OBJDIR) && rm * || true

install: $(EXECUTABLE)
	install -m 0755 $(EXECUTABLE) $(CODEROOT)/shared/bin

remove:
	rm  $(CODEROOT)/shared/bin/$(EXECUTABLE)

test: $(TOBJECT)
	$(CC) $(TOBJECT) $(LDFLAGS) -o $(TSTDIR)/$@
