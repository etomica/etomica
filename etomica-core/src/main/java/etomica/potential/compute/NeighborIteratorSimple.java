package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.BondingInfo;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

public class NeighborIteratorSimple implements NeighborIterator {

    private final Box box;
    private final BondingInfo bondingInfo;
    private final boolean isPureAtoms;
    private final Space space;

    public NeighborIteratorSimple(Box box, BondingInfo bondingInfo, boolean isPureAtoms) {
        this.box = box;
        this.bondingInfo = bondingInfo;
        this.isPureAtoms = isPureAtoms;
        this.space = box.getSpace();
    }

    @Override
    public void iterUpNeighbors(int iAtom, NeighborConsumer consumer) {
        IAtomList atoms = box.getLeafList();
        Boundary boundary = box.getBoundary();

        for (int i = iAtom + 1; i < atoms.size(); i++) {
            IAtom atom1 = atoms.get(iAtom);
            IAtom atom2 = atoms.get(i);

            if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                continue;
            }

            Vector dr = space.makeVector();
            dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
            boundary.nearestImage(dr);
            consumer.accept(atom2, dr);
        }
    }

    @Override
    public void iterDownNeighbors(int iAtom, NeighborConsumer consumer) {
        IAtomList atoms = box.getLeafList();
        Boundary boundary = box.getBoundary();

        for (int i = 1; i < iAtom; i++) {
            IAtom atom1 = atoms.get(iAtom);
            IAtom atom2 = atoms.get(i);

            if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                continue;
            }

            Vector dr = space.makeVector();
            dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
            boundary.nearestImage(dr);
            consumer.accept(atom2, dr);
        }
    }

    @Override
    public void iterAllNeighbors(int iAtom, NeighborConsumer consumer) {
        IAtomList atoms = box.getLeafList();
        Boundary boundary = box.getBoundary();

        for (int i = 0; i < atoms.size(); i++) {
            if (i == iAtom) { continue; }
            IAtom atom1 = atoms.get(iAtom);
            IAtom atom2 = atoms.get(i);

            if (bondingInfo.skipBondedPair(isPureAtoms, atom1, atom2)) {
                continue;
            }

            Vector dr = space.makeVector();
            dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
            boundary.nearestImage(dr);
            consumer.accept(atom2, dr);
        }
    }
}
