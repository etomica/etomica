/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell.molecule;

import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.MoleculeSetSinglet;
import etomica.atom.iterator.MoleculeIterator;
import etomica.atom.iterator.MoleculeIteratorArrayListSimple;
import etomica.box.BoxAgentManager;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLattice;

/**
 * Returns occupants of all cells as iterates.

 * @author Tai Boon Tan
 *
 */
public class MoleculeIteratorCell implements MoleculeIterator, java.io.Serializable {

    /**
     * Constructor makes iterator that must have box specified and then be
     * reset() before iteration.
     * 
     * @param D
     *            the dimension of the space of the simulation (used to
     *            construct cell iterators)
     */
	public MoleculeIteratorCell(int D, BoxAgentManager agentManager) {
        cellIterator = new RectangularLattice.Iterator(D);
        moleculeIterator = new MoleculeIteratorArrayListSimple();
        boxAgentManager = agentManager;
        moleculeSetSinglet = new MoleculeSetSinglet();
	}

	public void setBox(Box box) {
        CellLattice lattice = ((NeighborCellManagerMolecular)boxAgentManager.getAgent(box)).getLattice();
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
    
    public final IMoleculeList next() {
        moleculeSetSinglet.atom = nextMolecule();
        return moleculeSetSinglet;
    }
    
    public IMolecule nextMolecule() {
        IMolecule nextMolecule = moleculeIterator.nextMolecule();
        while (nextMolecule == null) {
            if(cellIterator.hasNext()) {
                moleculeIterator.setList(((CellMolecular)cellIterator.next()).occupants());
                moleculeIterator.reset();
            } else {//no more cells at all
                break;
            }
            nextMolecule = moleculeIterator.nextMolecule();
        }
        return nextMolecule;
    }
    
    public void unset() {
        moleculeIterator.unset();
    }

    /**
     * Returns 1, indicating that this is an atom iterator.
     */
    public int nBody() {
        return 1;
    }
    
    public void reset() {
        cellIterator.reset();
        moleculeIterator.unset();
    }
    
    private static final long serialVersionUID = 1L;
    private final MoleculeIteratorArrayListSimple moleculeIterator;
    private final RectangularLattice.Iterator cellIterator;
    private final BoxAgentManager boxAgentManager;
    protected final MoleculeSetSinglet moleculeSetSinglet;
}
