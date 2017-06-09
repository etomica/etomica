/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell.molecule;

import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.MoleculePair;
import etomica.atom.MoleculeToMoleculeSetFixed;
import etomica.atom.iterator.IteratorDirective;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.atom.iterator.MoleculeIterator;
import etomica.atom.iterator.MoleculeIteratorArrayList;
import etomica.atom.iterator.MoleculeIteratorArrayListSimple;
import etomica.atom.iterator.MoleculesetIteratorPDT;
import etomica.box.BoxAgentManager;
import etomica.lattice.CellLattice;

/**
 * Generates pairs that are cell-based neighbors of a specific Molecule. Iteration is
 * performed using cell lists, which defines the neighboring molecules.
 * Direction is related to ordering of cells and, within a cell, ordering of
 * molecules in cell's occupant list.
 * 
 * @author Tai Boon Tan
 *
 */
public class Mpi1ACell implements MoleculesetIteratorPDT, MoleculesetIteratorCellular,
                                  java.io.Serializable {
    
    /**
     * Constructor makes iterator that must have box specified and then be
     * reset() before iteration.
     * 
     * @param D
     *            the dimension of the space of the simulation (used to
     *            construct cell iterators)
     * @param range
     *            the distance within which pairs of atoms are considered
     *            neighbors. Used to define neighbor cells; some iterates may
     *            exceed this separation
     *  
     */
    public Mpi1ACell(int D, double range, BoxAgentManager agentManager) {
        neighborIterator = new CellLattice.NeighborIterator(D, range);
        aiSeq = new MoleculeIteratorArrayListSimple();
        //this iterator is used to loop through list of occupants of atoms's cell;
        //construct with AtomToLinker that gives appropriate linker
//        MyAtomToLinker atomToLinker = new MyAtomToLinker();
       
        molToMolSetFixed = new MoleculeToMoleculeSetFixed();
        aiSeqDirectableUp = new MoleculeIteratorArrayList(IteratorDirective.Direction.UP, 1, molToMolSetFixed, molToMolSetFixed);
        aiSeqDirectableDn = new MoleculeIteratorArrayList(IteratorDirective.Direction.DOWN, 1, molToMolSetFixed, molToMolSetFixed);
        latticeIndex = new int[D];
        periodicity = new boolean[D];

        neighborIterator.setDirection(null);
        boxAgentManager = agentManager;
	}

	public void setBox(Box box) {
        cellManager = (NeighborCellManagerMolecular)boxAgentManager.getAgent(box);
        lattice = cellManager.getLattice();
        neighborIterator.setLattice(lattice);
        for (int i=0; i<periodicity.length; i++) {
            periodicity[i] = box.getBoundary().getPeriodicity(i);
        }
        neighborIterator.setPeriodicity(periodicity);
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

    public IMoleculeList next() {
        IMolecule innerMolecule = aiInner.nextMolecule();
        if (innerMolecule == null) {
            innerMolecule = advanceLists();
            if (innerMolecule == null) {
                return null;
            }
        }
        if (upListNow) {
            pair.atom1 = innerMolecule;
            pair.atom0 = targetMolecule;
        }
        else {
            pair.atom0 = innerMolecule;
            pair.atom1 = targetMolecule;
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
        if(targetMolecule == null) {
            unset();
            return;
        }
        inCentralCell = true;
        upListNow = (direction != IteratorDirective.Direction.DOWN);
        doGoDown = (direction != IteratorDirective.Direction.UP);
        neighborIterator.checkDimensions();
        CellMolecular centralCell = cellManager.getCell(targetMolecule);
        lattice.latticeIndex(centralCell.latticeArrayIndex,latticeIndex);
        neighborIterator.setSite(latticeIndex);
        neighborIterator.setDirection(upListNow ? IteratorDirective.Direction.UP : IteratorDirective.Direction.DOWN);
        neighborIterator.reset();
        
        //start with targetMolecule's cell
        molToMolSetFixed.setArrayList(centralCell.occupants());
        if (upListNow) {
            aiSeqDirectableUp.setMolecule(targetMolecule);
            aiSeqDirectableUp.reset();
            pair.atom0 = targetMolecule;
            aiInner = aiSeqDirectableUp;
            return;
        }
        aiSeqDirectableDn.setMolecule(targetMolecule);
        aiSeqDirectableDn.reset();
        pair.atom1 = targetMolecule;
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
    public void setTarget(IMolecule newTargetMolecule) {
        targetMolecule = newTargetMolecule;
    }

    // Moves to next neighbor-cell list that can provide an iterate
    // This should be invoked only if aiInner.hasNext is false
    private IMolecule advanceLists() {
        if (inCentralCell && upListNow && doGoDown) {
            aiSeqDirectableDn.setMolecule(targetMolecule);
            aiSeqDirectableDn.reset();
            upListNow = false;
            aiInner = aiSeqDirectableDn;
            IMolecule molecule = aiSeqDirectableDn.nextMolecule();
            if (molecule != null) {
                return molecule;
            }
        }
        if (direction == null && inCentralCell) {
            upListNow = true;
        }
        inCentralCell = false;
        aiInner = aiSeq;//need to switch from aiSeqDirectableXX
        do {
            //advance neighbor cell 
            if(neighborIterator.hasNext()) {
                aiSeq.setList(((CellMolecular)neighborIterator.next()).occupants());
                aiSeq.reset();
            } else if (upListNow && doGoDown) {
                // handle "down" cells
                upListNow = false;
                neighborIterator.setDirection(IteratorDirective.Direction.DOWN);
                neighborIterator.reset();
                aiSeq.setList(((CellMolecular)neighborIterator.next()).occupants());
                aiSeq.reset();
            } else {
                //no more cells
                return null;
            }
            IMolecule molecule = aiSeq.nextMolecule();
            if (molecule != null) {
                return molecule;
            }
        } while(true);
    }

    /**
     * @return Returns the cellIterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator() {
        return neighborIterator;
    }
   
    private static final long serialVersionUID = 1L;
    private final CellLattice.NeighborIterator neighborIterator;
    private final MoleculeIteratorArrayList aiSeqDirectableUp, aiSeqDirectableDn;
    private final MoleculeIteratorArrayListSimple aiSeq;
    private final MoleculeToMoleculeSetFixed molToMolSetFixed;
    private final MoleculePair pair = new MoleculePair();
    private final int[] latticeIndex;
    private IteratorDirective.Direction direction;
    private boolean doGoDown, upListNow;
    private boolean inCentralCell;
    private IMolecule targetMolecule;
    private final BoxAgentManager boxAgentManager;
    private NeighborCellManagerMolecular cellManager;
    protected final boolean[] periodicity;
    
    private CellLattice lattice;
    
    private MoleculeIterator aiInner;

}
