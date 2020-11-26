package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.integrator.IntegratorListener;
import etomica.potential.IPotentialAtomic;

public interface NeighborManager {

    NeighborIterator makeNeighborIterator();

    default void setPairPotentials(IPotentialAtomic[][] potentials) {
    }

    void init();

    default void setPotentialRange(double range) {}

    IntegratorListener makeIntegratorListener();

    void updateAtom(IAtom atom);
}
