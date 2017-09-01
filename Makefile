ALL := libCalculus libLinkedList libFormaldehyde
ALL += HOCH_MolecularDynamics_MPI HOCH_MicroCanonical HOCH_Geodesics_MPI

all: libs execs

.PHONY: all execs libs clean remove $(ALL)

execs: HOCH_Geodesics_MPI HOCH_MicroCanonical HOCH_MolecularDynamics_MPI

libs: libFormaldehyde libLinkedList libCalculus

HOCH_Geodesics_MPI: libLinkedList libFormaldehyde
	$(MAKE) -C $@ install

HOCH_MicroCanonical: libLinkedList libFormaldehyde
	$(MAKE) -C $@ install

HOCH_MolecularDynamics_MPI: libLinkedList libFormaldehyde
	$(MAKE) -C $@ install

libFormaldehyde: libCalculus
	$(MAKE) -C $@ install

libLinkedList:
	$(MAKE) -C $@ install

libCalculus:
	$(MAKE) -C $@ install

clean:
	$(foreach dep, $(ALL), $(MAKE) -C $(dep) clean;)

remove:
	$(foreach dep, $(ALL), $(MAKE) -C $(dep) remove || true;)
