CODEROOT=..
SRCDIR=./src
OBJDIR=./build
TSTDIR=./test

LIBRARY := Calculus
TSOURCES := $(TSTDIR)/test.c

SHELL = /bin/sh
CC := gcc-6.2
#CC := mpixlc_r

CFLAGS := -g -O2 -DHAVE_INLINE
WARNINGS := -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Werror -Wuninitialized -Wconversion
CFLAGS += -I$(HOME)/toor/include -I$(CODEROOT)/shared/include
CFLAGS += -std=c99 -fshort-enums -fno-common
CFLAGS += $(WARNINGS)

LDFLAGS := -lgsl -lgslcblas -lblas -lm
LDFLAGS += -L$(HOME)/toor/lib
# LDFLAGS += -L$(HOME)/cross/lib -L$(CODEROOT)/shared/lib
# LDFLAGS += -L/soft/libraries/alcf/current/xl/CBLAS/lib
# LDFLAGS += -L/soft/libraries/alcf/current/xl/BLAS/lib
# LDFLAGS += -L/soft/apps/ibmcmp-feb2013/xlf/bg/14.1/bglib64 -lxlopt -lxl -lxlf90_r -lxlfmath
# LDFLAGS += -s

LIBNAME := lib$(LIBRARY).a

HEADERS := $(wildcard $(SRCDIR)/*.h)
SOURCES := $(wildcard $(SRCDIR)/*.c)
OBJECTS := $(patsubst $(SRCDIR)%.c,$(OBJDIR)%.o, $(SOURCES))

TOBJECTS := $(patsubst $(TSTDIR)%.c,$(OBJDIR)%.o, $(TSOURCES))

all: $(SOURCES) $(HEADERS) $(LIBRARY)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/%.o : $(TSTDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(LIBRARY): $(OBJECTS) $(HEADERS)
	ar rs $(LIBNAME) $(OBJECTS)

test: $(LIBRARY) $(TOBJECTS)
	$(CC) $(TOBJECTS) -L./ -l$(LIBRARY) $(LDFLAGS) -o $(TSTDIR)/$@
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

