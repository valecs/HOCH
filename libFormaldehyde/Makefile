CODEROOT=..
SRCDIR=./src
OBJDIR=./build
TSTDIR=./test
FSRCDIR=$(SRCDIR)/fortran

LIBRARY := Formaldehyde
SHELL = /bin/sh

FFLAGS := -O3
CFLAGS := -g -O2 -DHAVE_INLINE
CFLAGS += -I$(HOME)/toor/include -I$(HOME)/local/include -I$(CODEROOT)/shared/include

WARNINGS := -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Werror -Wuninitialized -Wconversion

LDFLAGS := -lgsl -lgslcblas -lblas -lm
LDFLAGS += -L$(CODEROOT)/shared/lib -L$(HOME)/toor/lib

XLKEY := /wrappers/xl
ifeq ($(XLKEY),$(findstring $(XLKEY), $(PATH)))
	CC := mpixlc_r
	FC := mpixlf77_r
#	since xl compilers don't mangle fortran names in the same was as gcc
	CFLAGS += -DEXTGETPOT=getpot
	CFLAGS += -qlanglvl=stdc99
	LDFLAGS += -L$(HOME)/cross/lib
	LDFLAGS += -L/soft/libraries/alcf/current/xl/CBLAS/lib
	LDFLAGS += -L/soft/libraries/alcf/current/xl/BLAS/lib
	LDFLAGS += -L/soft/apps/ibmcmp-feb2013/xlf/bg/14.1/bglib64 -lxlopt -lxl -lxlf90_r -lxlfmath
	LDFLAGS += -s
else # Assume GCC or compatible
	CC := mpicc
	FC := gfortran
	CFLAGS += -std=c99 -fshort-enums -fno-common
	CFLAGS += $(WARNINGS)
endif

LIBNAME := lib$(LIBRARY).a

HEADERS := $(wildcard $(SRCDIR)/*.h)
SOURCES := $(wildcard $(SRCDIR)/*.c) 
OBJECTS := $(patsubst $(SRCDIR)%.c,$(OBJDIR)%.o, $(SOURCES))

FSOURCES := $(FSRCDIR)/bonds.f $(FSRCDIR)/getpot.f $(FSRCDIR)/h2co-ccsdt.f $(FSRCDIR)/h2co-mrci.f
FOBJECTS := $(patsubst $(FSRCDIR)%.f,$(OBJDIR)%.o, $(FSOURCES))

TSOURCES := $(wildcard $(TSTDIR)/*.c) 
TOBJECTS := $(patsubst $(TSTDIR)%.c,$(OBJDIR)%.o, $(TSOURCES))

all: $(SOURCES) $(HEADERS) $(LIBRARY)

$(OBJDIR)/%.o : $(FSRCDIR)/%.f
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/%.o : $(TSTDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(LIBRARY): $(FOBJECTS) $(OBJECTS) $(HEADERS)
	ar rs $(LIBNAME) $(OBJECTS) $(FOBJECTS)

cobjects: $(OBJECTS)

test: $(LIBRARY) $(TOBJECTS)
	$(CC) $(TOBJECTS) -L./ -L$(CODEROOT)/shared/lib -l$(LIBRARY) -lCalculus $(LDFLAGS) -o $(TSTDIR)/$@
	./$(TSTDIR)/$@ > ./$(TSTDIR)/$@.out

clean:
	rm $(LIBNAME) || true
	cd $(TSTDIR) && rm test test.out || true
	cd $(OBJDIR) && rm * || true

install: $(LIBRARY) $(HEADERS)
	install -m 444 $(LIBNAME) $(CODEROOT)/shared/lib
	install -m 444 $(HEADERS) $(CODEROOT)/shared/include

remove:
	rm $(CODEROOT)/shared/lib/$(LIBNAME)
	$(foreach H,$(patsubst $(SRCDIR)%.h,$(CODEROOT)/shared/include%.h, $(HEADERS)),rm $(H);)

