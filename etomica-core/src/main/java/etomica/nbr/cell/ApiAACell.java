/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.cell;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.atom.AtomArrayList;
import etomica.atom.iterator.ApiInterArrayList;
import etomica.atom.iterator.ApiIntraArrayList;
import etomica.atom.iterator.AtomLeafsetIterator;
import etomica.atom.iterator.IteratorDirective;
import etomica.lattice.CellLattice;
import etomica.lattice.RectangularLattice;

/**
 * Returns iterates formed from all cell-based neighbor pairs.
 */
public class ApiAACell implements AtomsetIteratorCellular, java.io.Serializable {

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
     */
	public ApiAACell(int D, double range, Box box) {
        cellIterator = new RectangularLattice.Iterator(D);
        neighborIterator = new CellLattice.NeighborIterator(D, range);
        neighborIterator.setDirection(IteratorDirective.Direction.UP);
        interListIterator = new ApiInterArrayList(new AtomArrayList(), new AtomArrayList());
        intraListIterator = new ApiIntraArrayList();
        listIterator = intraListIterator;
        periodicity = new boolean[D];
        this.box = box;
	}

	public void setLattice(CellLattice lattice) {
        cellIterator.setLattice(lattice);
		neighborIterator.setLattice(lattice);
        unset();
	}
    
	/**
     * Returns the number of atom pairs the iterator will return if reset and
     * iterated in its present state.
     */
	public int size() {
        int count = 0;
        reset();
        for (Object a = nextPair(); a != null; a = nextPair()) {
            count++;
        }
        return count;
	}
	
    public IAtomList nextPair() {
        IAtomList nextPair = listIterator.next();
        if (nextPair == null) {
            return advanceLists();
        }
        return nextPair;
    }
    
    public void unset() {
        listIterator.unset();
    }

    /**
     * Returns 2, indicating that this is a pair iterator.
     */
    public int nBody() {
        return 2;
    }
    
    public void reset() {
        for (int i=0; i<periodicity.length; i++) {
            periodicity[i] = box.getBoundary().getPeriodicity(i);
        }
        neighborIterator.setPeriodicity(periodicity);
        cellIterator.reset();
        neighborIterator.checkDimensions();
        neighborIterator.unset();
        listIterator.unset();
        //System.out.println("reset in ApiAACell");
    }//end of reset
    
    // Moves to next pair of lists that can provide an iterate
    // This should be invoked only if listIterator.hasNext is false
    private IAtomList advanceLists() {
        do {
              //advance neighbor cell
            if(neighborIterator.hasNext()) {
                interListIterator.setInnerList(((Cell)neighborIterator.next()).occupants());
                listIterator = interListIterator;
                interListIterator.reset();
                IAtomList pair = listIterator.next();
                if (pair != null) {
                    return pair;
                }

                //advance central cell and set up neighbor cell iterator if
                // central cell has some molecules
            } else if(cellIterator.hasNext()) {
                AtomArrayList list = ((Cell)cellIterator.peek()).occupants();
                neighborIterator.setSite(cellIterator.nextIndex());

                if(!list.isEmpty()) {//central cell has molecules
                    interListIterator.setOuterList(list); //for neighbor-cell looping
                    intraListIterator.setList(list);//for intra-cell looping
                    neighborIterator.reset();

                    listIterator = intraListIterator;
                    intraListIterator.reset();
                        
                    IAtomList pair = listIterator.next();
                    if (pair != null) {
                        return pair;
                    }
                } else {//no molecules in central cell
                    neighborIterator.unset();
                    listIterator.unset();
                }
            } else {//no more cells at all
                return null;
            }
        } while(true);
    }//end of advanceCell
    
    /**
     * @return Returns the cellIterator.
     */
    public CellLattice.NeighborIterator getNbrCellIterator() {
        return neighborIterator;
    }
   
    private static final long serialVersionUID = 1L;
    private AtomLeafsetIterator listIterator;
    private final Box box;
    private final ApiIntraArrayList intraListIterator;
    private final ApiInterArrayList interListIterator;
    private final CellLattice.NeighborIterator neighborIterator;
    private final RectangularLattice.Iterator cellIterator;
    protected final boolean[] periodicity;
}
