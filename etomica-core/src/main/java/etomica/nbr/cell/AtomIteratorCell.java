/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.AtomSetSinglet;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.BoxAgentManager;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLattice;

/**
 * Returns occupants of all cells as iterates.
 */

public class AtomIteratorCell implements AtomIterator, java.io.Serializable {

    /**
     * Constructor makes iterator that must have box specified and then be
     * reset() before iteration.
     * 
     * @param D
     *            the dimension of the space of the simulation (used to
     *            construct cell iterators)
     */
	public AtomIteratorCell(int D, BoxAgentManager agentManager) {
        cellIterator = new RectangularLattice.Iterator(D);
        atomIterator = new AtomIteratorArrayListSimple();
        boxAgentManager = agentManager;
        atomSetSinglet = new AtomSetSinglet();
	}

	public void setBox(Box box) {
        CellLattice lattice = ((NeighborCellManager)boxAgentManager.getAgent(box)).getLattice();
        cellIterator.setLattice(lattice);
        unset();
	}
    
	/**
     * Returns the number of atoms the iterator will return if reset and
     * iterated in its present state.
     */
	public int size() {
        int count = 0;
        reset();
        for (Object a = next(); a != null; a = next()) {
            count++;
        }
        return count;
	}
	
    public boolean hasNext() {
        throw new RuntimeException("jfdka");
    }
    
    public final IAtomList next() {
        atomSetSinglet.atom = nextAtom();
        return atomSetSinglet;
    }
    
    public IAtom nextAtom() {
        IAtom nextAtom = atomIterator.nextAtom();
        while (nextAtom == null) {
            if(cellIterator.hasNext()) {
                atomIterator.setList(((Cell)cellIterator.next()).occupants());
                atomIterator.reset();
            } else {//no more cells at all
                break;
            }
            nextAtom = atomIterator.nextAtom();
        }
        return nextAtom;
    }
    
    public void unset() {
        atomIterator.unset();
    }

    /**
     * Returns 1, indicating that this is an atom iterator.
     */
    public int nBody() {
        return 1;
    }
    
    public void reset() {
        cellIterator.reset();
        atomIterator.unset();
    }
    
    private static final long serialVersionUID = 1L;
    private final AtomIteratorArrayListSimple atomIterator;
    private final RectangularLattice.Iterator cellIterator;
    private final BoxAgentManager boxAgentManager;
    protected final AtomSetSinglet atomSetSinglet;
}
