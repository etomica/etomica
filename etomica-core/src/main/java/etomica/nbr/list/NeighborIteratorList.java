package etomica.nbr.list;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.NeighborIterator;
import etomica.space.Space;
import etomica.space.Vector;

public class NeighborIteratorList implements NeighborIterator {
    private final NeighborListManagerFasterer nbrManager;
    private final Box box;
    private final Space space;

    public NeighborIteratorList(NeighborListManagerFasterer nbrManager, Box box) {
        this.nbrManager = nbrManager;
        this.box = box;
        this.space = box.getSpace();
    }

    @Override
    public void iterUpNeighbors(int iAtom, NeighborConsumer consumer) {
        IAtomList atoms = box.getLeafList();
        IAtom atom1 = atoms.get(iAtom);
        Vector ri = atom1.getPosition();
        int iNumNbrs = nbrManager.numAtomNbrsUp[iAtom];
        int[] iNbrs = nbrManager.nbrs[iAtom];
        Vector[] iNbrBoxOffsets = nbrManager.nbrBoxOffsets[iAtom];

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
    }

    @Override
    public void iterDownNeighbors(int iAtom, NeighborConsumer consumer) {
        IAtomList atoms = box.getLeafList();
        IAtom atom1 = atoms.get(iAtom);
        Vector ri = atom1.getPosition();
        int iNumNbrs = nbrManager.numAtomNbrsDn[iAtom];
        int[] iNbrs = nbrManager.nbrs[iAtom];
        Vector[] iNbrBoxOffsets = nbrManager.nbrBoxOffsets[iAtom];
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

    @Override
    public void iterAllNeighbors(int iAtom, NeighborConsumer consumer) {
        throw new UnsupportedOperationException();
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
