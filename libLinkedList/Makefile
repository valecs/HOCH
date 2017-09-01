CODEROOT= ..
SRCDIR=./src
OBJDIR=./build
TSTDIR=./test

LIBRARY := LinkedList

SHELL = /bin/sh

CFLAGS :=  -g -O2
WARNINGS := -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Werror -Wuninitialized -Wconversion

XLKEY := /wrappers/xl
ifeq ($(XLKEY),$(findstring $(XLKEY), $(PATH)))
	CC := mpixlc_r
	CFLAGS += -qlanglvl=stdc99
else # Assume GCC or compatible
	CC := gcc
	CFLAGS += -std=c99 -fshort-enums -fno-common
	CFLAGS += $(WARNINGS)
endif

SRCDIR:= ./src
OBJDIR:= ./build
TSTDIR:= ./test

LIBNAME := lib$(LIBRARY).a
HEADERS := $(wildcard $(SRCDIR)/*.h)
SOURCES := $(wildcard $(SRCDIR)/*.c)
OBJECTS := $(patsubst $(SRCDIR)%.c,$(OBJDIR)%.o, $(SOURCES))

TSOURCES := $(wildcard $(TSTDIR)/*.c) 
TOBJECTS := $(patsubst $(TSTDIR)%.c,$(OBJDIR)%.o, $(TSOURCES))

all: $(SOURCES) $(HEADERS) $(LIBRARY)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/%.o : $(TSTDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

$(LIBRARY): $(OBJECTS) $(HEADERS)
	ar rs $(LIBNAME) $(OBJECTS)

test: $(LIBRARY) $(TOBJECTS)
	$(CC) $(TOBJECTS) $(LDFLAGS) -L./ -l$(LIBRARY) -o $(TSTDIR)/$@
	./$(TSTDIR)/$@ > ./$(TSTDIR)/$@.out

clean:
	rm $(LIBNAME) || true
	cd $(TSTDIR) && rm test || true
	cd $(OBJDIR) && rm * || true

install: $(LIBRARY) $(HEADERS)
	install -m 444 $(LIBNAME) $(CODEROOT)/shared/lib
	install -m 444 $(HEADERS) $(CODEROOT)/shared/include

remove:
	rm $(CODEROOT)/shared/lib/$(LIBNAME)
	$(foreach H,$(patsubst $(SRCDIR)%.h,$(CODEROOT)/shared/include%.h, $(HEADERS)),rm $(H);)

