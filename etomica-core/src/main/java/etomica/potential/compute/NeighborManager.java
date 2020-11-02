package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.integrator.IntegratorListener;
import etomica.potential.Potential2Soft;

public interface NeighborManager {

    NeighborIterator makeNeighborIterator();

    default void setPairPotentials(Potential2Soft[][] potentials) {}

    void init();

    default void setPotentialRange(double range) {}

    IntegratorListener makeIntegratorListener();

    void updateAtom(IAtom atom);
}
