package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.NeighborConsumerHard;
import etomica.potential.compute.NeighborIteratorHard;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Neighbor iterator that uses neighbor lists and can pass the pair state to
 * the callback for hard MD.
 */
public class NeighborIteratorListHard implements NeighborIteratorHard {
    private final NeighborListManagerFastererHard nbrManager;
    private final Box box;
    private final Space space;

    public NeighborIteratorListHard(NeighborListManagerFastererHard nbrManager, Box box) {
        this.nbrManager = nbrManager;
        this.box = box;
        this.space = box.getSpace();
    }

    @Override
    public void iterUpNeighbors(int iAtom, NeighborConsumerHard consumer, double falseTime) {
        iterUpNeighbors(iAtom, consumer);
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
    public void iterDownNeighbors(int iAtom, NeighborConsumerHard consumer, double falseTime) {
        iterDownNeighbors(iAtom, consumer);
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
    public void iterAllNeighbors(int iAtom, NeighborConsumerHard consumer, double falseTime) {
        iterAllNeighbors(iAtom, consumer);
    }

    @Override
    public void iterAllNeighbors(int iAtom, NeighborConsumer consumer) {
        IAtomList atoms = box.getLeafList();
        IAtom atom1 = atoms.get(iAtom);
        Vector ri = atom1.getPosition();
        int iNumNbrs = nbrManager.numAtomNbrsUp[iAtom];
        int[] iNbrs = nbrManager.nbrs[iAtom];
        Vector[] iNbrBoxOffsets = nbrManager.nbrBoxOffsets[iAtom];
        double sum = 0;

        for (int j = 0; j < iNumNbrs; j++) {
            int jAtom = iNbrs[j];
            IAtom atom2 = atoms.get(jAtom);
            Vector rj = atom2.getPosition();
            Vector jbo = iNbrBoxOffsets[j];
            Vector rij = space.makeVector();
            rij.Ev1Mv2(rj, ri);
            rij.PE(jbo);
            consumer.accept(atom2, rij);
        }

        iNumNbrs = nbrManager.numAtomNbrsDn[iAtom];
        iNbrs = nbrManager.nbrs[iAtom];
        int maxNbrs = iNbrs.length;

        for (int j = maxNbrs - 1; j > maxNbrs - 1 - iNumNbrs; j--) {
            int jAtom = iNbrs[j];
            IAtom atom2 = atoms.get(jAtom);
            Vector rj = atom2.getPosition();
            Vector jbo = iNbrBoxOffsets[j];
            Vector rij = space.makeVector();
            rij.Ev1Mv2(rj, ri);
            rij.ME(jbo);
            consumer.accept(atom2, rij);
        }
    }

    public double iterAndSumAllNeighbors(IAtom atom1, SuperNbrConsumer consumer) {
        IAtomList atoms = box.getLeafList();
        int iAtom = atom1.getLeafIndex();
        Vector ri = atom1.getPosition();
        int iNumNbrs = nbrManager.numAtomNbrsUp[iAtom];
        int[] iNbrs = nbrManager.nbrs[iAtom];
        Vector[] iNbrBoxOffsets = nbrManager.nbrBoxOffsets[iAtom];
        double sum = 0;

        for (int j = 0; j < iNumNbrs; j++) {
            int jAtom = iNbrs[j];
            IAtom atom2 = atoms.get(jAtom);
            Vector rj = atom2.getPosition();
            Vector jbo = iNbrBoxOffsets[j];
            Vector rij = space.makeVector();
            rij.Ev1Mv2(rj, ri);
            rij.PE(jbo);
            sum += consumer.accept(atom1, atom2, rij);
        }

        iNumNbrs = nbrManager.numAtomNbrsDn[iAtom];
        iNbrs = nbrManager.nbrs[iAtom];
        int maxNbrs = iNbrs.length;

        for (int j = maxNbrs - 1; j > maxNbrs - 1 - iNumNbrs; j--) {
            int jAtom = iNbrs[j];
            IAtom atom2 = atoms.get(jAtom);
            Vector rj = atom2.getPosition();
            Vector jbo = iNbrBoxOffsets[j];
            Vector rij = space.makeVector();
            rij.Ev1Mv2(rj, ri);
            rij.ME(jbo);
            sum += consumer.accept(atom1, atom2, rij);
        }

        return sum;
    }
}
