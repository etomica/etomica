/*
 * History
 * Created on Nov 23, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.AtomList;
import etomica.lattice.AbstractLattice;
import etomica.utility.Arrays;
import etomica.lattice.SiteFactory;

/**
 * Site used to form array of cells for cell-based neighbor listing.  Each
 * cell is capable of holding lists of atoms that are in them.
 */
public class NeighborCell {

    /**
     * Adds a new AtomList at the end of the array of occupant lists.
     */
    public void addOccupantList() {
        occupants = (AtomList[])Arrays.addObject(occupants,new AtomList());
    }
    /**
     * Removes the AtomList from and resizes the occupant array. Throws
     * ArrayIndexOutOfBoundsException if index exceeds largest array element.
     */
    public void removeOccupantList(int index) {
        occupants = (AtomList[])Arrays.removeObject(occupants,occupants[index]);
    }
    
    public AtomList[] occupants() {return occupants;}
    
    private AtomList[] occupants = new AtomList[0];

    public static final SiteFactory FACTORY = new SiteFactory() {
        public Object makeSite(AbstractLattice lattice, int[] coord) {
            return new NeighborCell();
        }
    };
}
