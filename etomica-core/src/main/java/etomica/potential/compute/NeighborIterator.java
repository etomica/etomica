package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.space.Vector;

public interface NeighborIterator {

    void iterUpNeighbors(int iAtom, NeighborConsumer consumer);

    void iterDownNeighbors(int iAtom, NeighborConsumer consumer);

    void iterAllNeighbors(int iAtom, NeighborConsumer consumer);

    default double iterAndSumAllNeighbors(IAtom atom1, SuperNbrConsumer consumer) {
        return 0;
    }

    /**
     * Interface for neighbor iteration callback.
     */
    interface NeighborConsumer {
        void accept(IAtom jAtom, Vector rij, int n);
    }

    interface SuperNbrConsumer {
        double accept(IAtom atom1, IAtom atom2, Vector rij, int n);
    }
}
