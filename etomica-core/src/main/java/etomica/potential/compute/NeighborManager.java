package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorListener;
import etomica.nbr.list.INeighborListener;
import etomica.potential.BondingInfo;
import etomica.potential.IPotential2;

public interface NeighborManager {

    NeighborIterator makeNeighborIterator();

    BondingInfo getBondingInfo();

    default void setPairPotentials(IPotential2[][] potentials) {
    }

    void init();

    default void setPotentialRange(double range) {
    }

    IntegratorListener makeIntegratorListener();

    void updateAtom(IAtom atom);

    Box getBox();

    interface NeighborEventSource {
        void addListener(INeighborListener listener);
    }
}
