/*
 * History
 * Created on Nov 19, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.Atom;
import etomica.AtomIteratorAtomDependent;
import etomica.AtomIteratorList;
import etomica.AtomLinker;
import etomica.AtomList;
import etomica.Species;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.lattice.NeighborManager;

/**
 * Iterator giving those atoms found in the neighboring cells of a
 * given atom.  Written specifically for case in which all relevant
 * atoms are at the molecule level of the atom hierarchy.
 */
//TODO make this AtomIteratorDirectable by filtering based on up & down
public class AtomIteratorNbrCell implements AtomIteratorAtomDependent {

    /**
     * Constructs iterator such that it will return molecules of the given
     * species as its iterates.
     */
    public AtomIteratorNbrCell(Species species, boolean upOnly) {
        nbrListIndex = species.getIndex();
        this.upOnly = upOnly;
    }

    /**
     * Specifies the atom whose neighbors will be returned as iterates.
     */
    public void setAtom(Atom atom) {
        this.atom = atom;
        NeighborCell referenceCell = ((Cellular)atom.seq).getCell();
        NeighborManager neighborManager = referenceCell.neighborManager();
        if (upOnly) {
            cellLinker = neighborManager.tab;
        }
        cellIterator.setList(neighborManager.neighbors());
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
        if (upOnly) {
            cellIterator.setFirst(cellLinker);
            cellIterator.reset();//reset cell iterator
            NeighborCell nextCell = (NeighborCell)cellIterator.nextAtom();
            AtomList[] occupants = nextCell.occupants();
            AtomList nbrList = occupants[nbrListIndex];
//            AtomList nbrList = ((NeighborCell)cellIterator.nextAtom()).occupants()[nbrListIndex];
            atomIterator.setList(nbrList);
            atomIterator.setFirst(((AtomSequencerCell)atom.seq).nbrLink);
            atomIterator.reset();
            atomIterator.nextAtom();
            if (!atomIterator.hasNext()) {
                advanceCell();
            }
        }
        else {
            cellIterator.reset();
            atomIterator.unset();
            advanceCell();//reset list iterator
        }
    }

    public void unset() {
        cellIterator.unset();
        atomIterator.unset();
    }

    // Moves to next cell that has an iterate
    private void advanceCell() {
        do {
            if(cellIterator.hasNext()) {
                AtomList nbrList = ((NeighborCell)cellIterator.nextAtom()).occupants()[nbrListIndex];
                atomIterator.setList(nbrList);
                atomIterator.reset();
                if (!upOnly && atomIterator.hasNext() && atomIterator.peek()[0] == atom) {
                    atomIterator.nextAtom();
                }
            } else {//no more cells at all
                break;
            }
        } while(!atomIterator.hasNext());
    }


    public Atom nextAtom() {
        Atom atom = atomIterator.nextAtom();
        if(!upOnly && atomIterator.peek()[0] == atom) atomIterator.next();//skip central atom
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
        if (upOnly) {
            cellIterator.setFirst(cellLinker);
        }
        cellIterator.reset();
        boolean firstCell = true;
        while(cellIterator.hasNext()) {
            AtomList nbrList = ((NeighborCell)cellIterator.nextAtom()).occupants()[nbrListIndex];
            atomIterator.setList(nbrList);
            if (upOnly && firstCell) {
                atomIterator.setFirst(((AtomSequencerCell)atom.seq).nbrLink);
                firstCell = false;
            }
            atomIterator.reset();
            while(atomIterator.hasNext()) {
                Atom[] next = atomIterator.next();
                if(!upOnly && next[0] != atom) action.actionPerformed(next);
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
    private final AtomIteratorList atomIterator = new AtomIteratorList();
    private final AtomIteratorList cellIterator = new AtomIteratorList();
    private final boolean upOnly;
    private AtomLinker cellLinker;
    
}
