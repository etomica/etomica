/*
 * History
 * Created on Nov 19, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.Atom;
import etomica.AtomIteratorAtomDependent;
import etomica.AtomIteratorListSimple;
import etomica.AtomList;
import etomica.Species;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.nbr.cell.IteratorFactoryCell.CellSequencer;

/**
 * Iterator giving those atoms found in the neighboring cells of a
 * given atom.  Written specifically for case in which all relevant
 * atoms are at the molecule level of the atom hierarchy.
 */

//TODO worry about direction issues

public class AtomIteratorNbrCell implements AtomIteratorAtomDependent {

    /**
     * Constructs iterator such that it will return molecules of the given
     * species as its iterates.
     */
    public AtomIteratorNbrCell(Species species) {
        nbrListIndex = species.getIndex();
    }

    /**
     * Specifies the atom whose neighbors will be returned as iterates.
     * Throws NullPointerException if null atom is specified.
     */
    public void setAtom(Atom atom) {
        this.atom = atom;
        NeighborCell referenceCell = ((Cellular)atom.seq).getCell();
        cellIterator.setList(referenceCell.neighborManager().neighbors());
    }

    public boolean contains(Atom[] atom) {
        if(atom == null || atom.length == 0) return false;

        AtomsetDetect detector = new AtomsetDetect(atom[0]);
        allAtoms(detector);
        return detector.detectedAtom();
    }

    public boolean hasNext() {
        return atomIterator.hasNext();
    }

    public void reset() {
        cellIterator.reset();//reset cell iterator
        advanceCell();//reset list iterator
    }

    public void unset() {
        cellIterator.unset();
        atomIterator.unset();
    }

    // Moves to next cell that has an iterate
    private void advanceCell() {
        do {
            if(cellIterator.hasNext()) {//cell iterator
                AtomList nbrList = ((NeighborCell)cellIterator.nextAtom()).occupants()[nbrListIndex];
                atomIterator.setList(nbrList);
                atomIterator.reset();
            } else {//no more cells at all
                break;
            }
        } while(!atomIterator.hasNext() || atomIterator.peek()[0] == atom);
    }


    public Atom nextAtom() {
        Atom atom = atomIterator.nextAtom();
        if(atomIterator.peek()[0] == atom) atomIterator.next();//skip central atom
        if(!atomIterator.hasNext()) advanceCell();
        return atom;
    }

    public Atom[] next() {
        atoms[0] = nextAtom();
        return atoms;
    }

    public Atom[] peek() {
        return atomIterator.peek();
    }

    public void allAtoms(AtomsetAction action) {
        if(atom == null) return;//atom will be null if setAtom not yet called
        cellIterator.reset();
        while(cellIterator.hasNext()) {
            AtomList nbrList = ((NeighborCell)cellIterator.nextAtom()).occupants()[nbrListIndex];
            atomIterator.setList(nbrList);
            atomIterator.reset();
            while(atomIterator.hasNext()) {
                Atom[] next = atomIterator.next();
                if(next[0] != atom) action.actionPerformed(next);
            }
        }
    }

    public int size() {
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
    }

    /**
     * Returns 1, indicating this is a single-atom iterator.
     */
    public int nBody() {
        return 1;
    }

    private Atom atom;
    private final int nbrListIndex;
    private final Atom[] atoms = new Atom[1];
    private final AtomIteratorListSimple atomIterator = new AtomIteratorListSimple();
    private final AtomIteratorListSimple cellIterator = new AtomIteratorListSimple();

}
