EXECUTABLE := HOCH_MicroCanonical

CODEROOT := ..
SRCDIR:= ./src
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

points: driver.sh $(EXECUTABLE)
	./driver.sh
