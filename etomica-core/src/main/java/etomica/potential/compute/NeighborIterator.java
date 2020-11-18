package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.nbr.cell.NeighborIteratorCellFasterer;
import etomica.space.Vector;

public interface NeighborIterator {

    void iterUpNeighbors(int iAtom, NeighborConsumer consumer);

    void iterDownNeighbors(int iAtom, NeighborConsumer consumer);

    void iterAllNeighbors(int iAtom, NeighborConsumer consumer);

    default double iterAndSumAllNeighbors(IAtom atom1, NeighborIteratorCellFasterer.SuperNbrConsumer consumer) {
        return 0;
    }

    interface NeighborConsumer {
        void accept(IAtom jAtom, Vector rij);
    }
}
