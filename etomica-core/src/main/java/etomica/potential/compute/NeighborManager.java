package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorListener;
import etomica.nbr.list.INeighborListener;
import etomica.potential.BondingInfo;
import etomica.potential.Potential2Soft;

public interface NeighborManager {

    NeighborIterator makeNeighborIterator();

    BondingInfo getBondingInfo();

    default void setPairPotentials(Potential2Soft[][] potentials) {
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
