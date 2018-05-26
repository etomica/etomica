/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell;

import etomica.atom.AtomPair;
import etomica.atom.AtomToAtomSetFixed;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.*;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.lattice.CellLattice;
import etomica.potential.IteratorDirective;
import etomica.potential.IteratorDirective.Direction;

/**
 * Generates pairs that are cell-based neighbors of a specific Atom. Iteration is
 * performed using cell lists, which defines the neighboring molecules.
 * Direction is related to ordering of cells and, within a cell, ordering of
 * molecules in cell's occupant list.
 */
public class Api1ACell implements AtomsetIteratorDirectable, AtomsetIteratorTargetable, AtomsetIteratorCellular {

    private final Box box;
    private final CellLattice.NeighborIterator neighborIterator;
    private final AtomIteratorArrayList aiSeqDirectableUp, aiSeqDirectableDn;
    private final AtomIteratorArrayListSimple aiSeq;
    private final AtomToAtomSetFixed atomToAtomSetFixed;
    private final AtomPair pair = new AtomPair();
    private final int[] latticeIndex;
    private IteratorDirective.Direction direction;
    private boolean doGoDown, upListNow;
    private boolean inCentralCell;
    private IAtom targetAtom;
    private final NeighborCellManager cellManager;
    private CellLattice lattice;
    private AtomIterator aiInner;
    /**
     * Constructor makes iterator that must have box specified and then be
     * reset() before iteration.
     *
     * @param D     the dimension of the space of the simulation (used to
     *              construct cell iterators)
     * @param range the distance within which pairs of atoms are considered
     *              neighbors. Used to define neighbor cells; some iterates may
     *              exceed this separation
     */
    public Api1ACell(double range, Box box, NeighborCellManager neighborCellManager) {
        this.box = box;
        cellManager = neighborCellManager;
        neighborIterator = new CellLattice.NeighborIterator(box.getSpace().D(), range);
        neighborIterator.setPeriodicity(box.getBoundary().getPeriodicity());
        lattice = cellManager.getLattice();
        neighborIterator.setLattice(lattice);
        aiSeq = new AtomIteratorArrayListSimple();
        //this iterator is used to loop through list of occupants of atoms's cell;
        //construct with AtomToLinker that gives appropriate linker
//        MyAtomToLinker atomToLinker = new MyAtomToLinker();
        atomToAtomSetFixed = new AtomToAtomSetFixed();
        aiSeqDirectableUp = new AtomIteratorArrayList(IteratorDirective.Direction.UP, 1, atomToAtomSetFixed, atomToAtomSetFixed);
        aiSeqDirectableDn = new AtomIteratorArrayList(IteratorDirective.Direction.DOWN, 1, atomToAtomSetFixed, atomToAtomSetFixed);
        latticeIndex = new int[box.getSpace().D()];

        neighborIterator.setDirection(null);
    }

    /**
     * Returns the number of atom pairs the iterator will return if
     * reset and iterated in its present state.
     */
    public int size() {
        int count = 0;
        reset();
        for (Object a = next(); a != null; a = next()) {
            count++;
        }
        return count;
    }

    public IAtomList next() {
        IAtom innerAtom = aiInner.nextAtom();
        if (innerAtom == null) {
            innerAtom = advanceLists();
            if (innerAtom == null) {
                return null;
            }
        }
        if (upListNow) {
            pair.atom1 = innerAtom;
            pair.atom0 = targetAtom;
        } else {
            pair.atom0 = innerAtom;
            pair.atom1 = targetAtom;
        }
        return pair;
    }

    public void unset() {
        aiInner.unset();
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public int nBody() {
        return 2;
    }

    public void reset() {
        if (targetAtom == null) {
            unset();
            return;
        }
        inCentralCell = true;
        upListNow = (direction != IteratorDirective.Direction.DOWN);
        doGoDown = (direction != IteratorDirective.Direction.UP);
        neighborIterator.checkDimensions();
        Cell centralCell = cellManager.getCell(targetAtom);
        lattice.latticeIndex(centralCell.latticeArrayIndex, latticeIndex);
        neighborIterator.setSite(latticeIndex);
        neighborIterator.setDirection(upListNow ? IteratorDirective.Direction.UP : IteratorDirective.Direction.DOWN);
        neighborIterator.reset();

        //start with targetMolecule's cell
        atomToAtomSetFixed.setArrayList(centralCell.occupants());
        if (upListNow) {
            aiSeqDirectableUp.setAtom(targetAtom);
            aiSeqDirectableUp.reset();
            pair.atom0 = targetAtom;
            aiInner = aiSeqDirectableUp;
            return;
        }
        aiSeqDirectableDn.setAtom(targetAtom);
        aiSeqDirectableDn.reset();
        pair.atom1 = targetAtom;
        aiInner = aiSeqDirectableDn;
    }

    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. Specification of a null direction indicates iteration in both directions
     * relative to the target. Direction is determined by ordering within occupant
     * list of cell of target atom, and then by the cell ordering of neighboring cells.
     */
    public void setDirection(Direction direction) {
        this.direction = direction;
        neighborIterator.setDirection(direction);
    }

    /**
     * Sets the target molecule with which all pairs are formed.  Molecule
     * is determined from the first atom of the array, which may be the molecule
     * itself or an atom that is part of it.  If the atom is null or is not
     * in one of the species given at construction, no iterates will be returned.
     */
    public void setTarget(IAtom newTargetAtom) {
        targetAtom = newTargetAtom;
    }

    // Moves to next neighbor-cell list that can provide an iterate
    // This should be invoked only if aiInner.hasNext is false
    private IAtom advanceLists() {
        if (inCentralCell && upListNow && doGoDown) {
            aiSeqDirectableDn.setAtom(targetAtom);
            aiSeqDirectableDn.reset();
            upListNow = false;
            aiInner = aiSeqDirectableDn;
            IAtom atom = aiSeqDirectableDn.nextAtom();
            if (atom != null) {
                return atom;
            }
        }
        if (direction == null && inCentralCell) {
            upListNow = true;
        }
        inCentralCell = false;
        aiInner = aiSeq;//need to switch from aiSeqDirectableXX
        do {
            //advance neighbor cell
            if (neighborIterator.hasNext()) {
                aiSeq.setList(((Cell) neighborIterator.next()).occupants());
                aiSeq.reset();
            } else if (upListNow && doGoDown) {
                // handle "down" cells
                upListNow = false;
                neighborIterator.setDirection(IteratorDirective.Direction.DOWN);
                neighborIterator.reset();
                aiSeq.setList(((Cell) neighborIterator.next()).occupants());
                aiSeq.reset();
            } else {
                //no more cells
                return null;
            }
            IAtom atom = aiSeq.nextAtom();
            if (atom != null) {
                return atom;
            }
        } while (true);
    }

    /**
     * @return Returns the cellIterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator() {
        return neighborIterator;
    }

}
