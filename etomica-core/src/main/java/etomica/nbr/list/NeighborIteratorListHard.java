package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.NeighborConsumerHard;
import etomica.potential.compute.NeighborIterator;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Neighbor iterator that uses neighbor lists and can pass the pair state to
 * the callback for hard MD.
 */
public class NeighborIteratorListHard implements NeighborIterator {
    private final NeighborListManagerFastererHard nbrManager;
    private final Box box;
    private final Space space;

    public NeighborIteratorListHard(NeighborListManagerFastererHard nbrManager, Box box) {
        this.nbrManager = nbrManager;
        this.box = box;
        this.space = box.getSpace();
    }

    @Override
    public void iterUpNeighbors(int iAtom, NeighborConsumer consumer) {
        final NeighborConsumerHard consumerHard = consumer instanceof NeighborConsumerHard ? (NeighborConsumerHard) consumer : null;
        IAtomList atoms = box.getLeafList();
        IAtom atom1 = atoms.get(iAtom);
        Vector ri = atom1.getPosition();
        int iNumNbrs = nbrManager.numAtomNbrsUp[iAtom];
        int[] iNbrs = nbrManager.nbrs[iAtom];
        Vector[] iNbrBoxOffsets = nbrManager.nbrBoxOffsets[iAtom];
        int[] iState = nbrManager.nbrState[iAtom];

        for (int j = 0; j < iNumNbrs; j++) {
            int jAtom = iNbrs[j];
            IAtom atom2 = atoms.get(jAtom);
            Vector rj = atom2.getPosition();
            Vector jbo = iNbrBoxOffsets[j];
            Vector rij = space.makeVector();
            rij.Ev1Mv2(rj, ri);
            rij.PE(jbo);
            if (consumerHard == null) {
                // we're like just computing energy
                consumer.accept(atom2, rij);
            } else {
                consumerHard.accept(atom2, rij, iState[j]);
            }
        }
    }

    @Override
    public void iterDownNeighbors(int iAtom, NeighborConsumer consumer) {
        final NeighborConsumerHard consumerHard = consumer instanceof NeighborConsumerHard ? (NeighborConsumerHard) consumer : null;
        IAtomList atoms = box.getLeafList();
        IAtom atom1 = atoms.get(iAtom);
        Vector ri = atom1.getPosition();
        int iNumNbrs = nbrManager.numAtomNbrsDn[iAtom];
        int[] iNbrs = nbrManager.nbrs[iAtom];
        Vector[] iNbrBoxOffsets = nbrManager.nbrBoxOffsets[iAtom];
        int maxNbrs = iNbrs.length;
        int[] iState = nbrManager.nbrState[iAtom];

        for (int j = maxNbrs - 1; j > maxNbrs - 1 - iNumNbrs; j--) {
            int jAtom = iNbrs[j];
            IAtom atom2 = atoms.get(jAtom);
            Vector rj = atom2.getPosition();
            Vector jbo = iNbrBoxOffsets[j];
            Vector rij = space.makeVector();
            rij.Ev1Mv2(rj, ri);
            rij.ME(jbo);
            if (consumerHard == null) {
                // we're like just computing energy
                consumer.accept(atom2, rij);
            } else {
                consumerHard.accept(atom2, rij, iState[j]);
            }
        }
    }

    @Override
    public void iterAllNeighbors(int iAtom, NeighborConsumer consumer) {
        throw new UnsupportedOperationException();
    }
}