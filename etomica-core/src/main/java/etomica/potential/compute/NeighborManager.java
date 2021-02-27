package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorListener;
import etomica.potential.BondingInfo;
import etomica.potential.IPotentialAtomic;

public interface NeighborManager {

    NeighborIterator makeNeighborIterator();

    BondingInfo getBondingInfo();

    default void setPairPotentials(IPotentialAtomic[][] potentials) {
    }

    void init();

    default void setPotentialRange(double range) {
    }

    IntegratorListener makeIntegratorListener();

    void updateAtom(IAtom atom);

    Box getBox();
}
