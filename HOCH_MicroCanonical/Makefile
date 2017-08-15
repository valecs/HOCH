EXECUTABLE := HOCH_MicroCanonical

CODEROOT := ..
SRCDIR:= ./src
OBJDIR:= ./build

SHELL = /bin/sh
CC := gcc
PGLINKER=pgcc

WARNINGS := -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Werror
WARNINGS += -Wconversion -Werror # gcc 4.1.2 still combines -Wconversion and -Wconversion-traditional
CFLAGS:=-std=c99 -fshort-enums -fno-common -g -DHAVE_INLINE -I/home/vale/toor/include -I$(CODEROOT)/shared/include -DDEBUG_HEAP
#CFLAGS += -DDEBUG
CFLAGS += $(WARNINGS)

LDFLAGS:=-lgsl -lm -lgslcblas -lgfortran -L/home/vale/toor/lib -L$(CODEROOT)/shared/lib
LDFLAGS += -lFormaldehyde-pgi -lCalculus -lLinkedList -lbstr
PGLDFLAGS = -L$/share/apps/pgi/linux86-64/11.4/lib -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl -lrt

HEADERS := $(wildcard $(SRCDIR)/*.h)
SOURCES := $(wildcard $(SRCDIR)/*.c)
OBJECTS := $(patsubst $(SRCDIR)%.c,$(OBJDIR)%.o, $(SOURCES))

all: $(SOURCES) $(HEADERS) $(EXECUTABLE)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) -O2 $< -o $@

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(PGLINKER) $(OBJECTS) $(LDFLAGS) $(PGLDFLAGS) -o $@
#	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

clean:
	rm $(EXECUTABLE) || true
	cd $(OBJDIR) && rm * || true

install: $(EXECUTABLE)
	install -m 0755 $(EXECUTABLE) $(CODEROOT)/shared/bin

remove:
	rm  $(CODEROOT)/shared/bin/$(EXECUTABLE)

points: driver.sh $(EXECUTABLE)
	./driver.sh
