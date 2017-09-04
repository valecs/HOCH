LIBS  := libCalculus libLinkedList libFormaldehyde
EXECS := HOCH_MolecularDynamics_MPI HOCH_MicroCanonical HOCH_Geodesics_MPI

all: libs execs

.PHONY: all execs libs clean remove

execs: $(EXECS)

libs: $(LIBS)

HOCH_Geodesics_MPI: shared/bin/HOCH_Geodesics_MPI

HOCH_MicroCanonical: shared/bin/HOCH_MicroCanonical

HOCH_MolecularDynamics_MPI: shared/bin/HOCH_MolecularDynamics_MPI

libFormaldehyde: shared/lib/libFormaldehyde.a

libLinkedList: shared/lib/libLinkedList.a

libCalculus: shared/lib/libCalculus.a

shared/lib/libCalculus.a:
	$(MAKE) -C $(basename $(notdir $@)) install

shared/lib/libFormaldehyde.a: libCalculus
	$(MAKE) -C $(basename $(notdir $@)) install

shared/lib/libLinkedList.a:
	$(MAKE) -C $(basename $(notdir $@)) install

shared/bin/HOCH_Geodesics_MPI: libLinkedList libFormaldehyde
	$(MAKE) -C $(notdir $@) install

shared/bin/HOCH_MicroCanonical: libLinkedList libFormaldehyde
	$(MAKE) -C $(notdir $@) install

shared/bin/HOCH_MolecularDynamics_MPI: libLinkedList libFormaldehyde
	$(MAKE) -C $(notdir $@) install

clean:
	$(foreach dep, $(LIBS),  $(MAKE) -C $(dep) clean;)
	$(foreach dep, $(EXECS), $(MAKE) -C $(dep) clean;)

remove:
	$(foreach dep, $(LIBS),  $(MAKE) -C $(dep) remove || true;)
	$(foreach dep, $(EXECS), $(MAKE) -C $(dep) remove || true;)
