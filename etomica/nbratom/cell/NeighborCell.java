/*
 * History
 * Created on Nov 23, 2004 by kofke
 */
package etomica.nbratom.cell;

import etomica.atom.AtomList;
import etomica.lattice.AbstractLattice;
import etomica.lattice.RectangularLattice;
import etomica.lattice.SiteFactory;

/**
 * Site used to form array of cells for cell-based neighbor listing.  Each
 * cell is capable of holding lists of atoms that are in them.
 */


public class NeighborCell {

    public NeighborCell(int latticeArrayIndex) {
        this.latticeArrayIndex = latticeArrayIndex;
    }
    
    public AtomList occupants() {return occupants;}
    
    private final AtomList occupants = new AtomList();
    final int latticeArrayIndex;//identifies site in lattice

    public static final SiteFactory FACTORY = new SiteFactory() {
        public Object makeSite(AbstractLattice lattice, int[] coord) {
            return new NeighborCell(((RectangularLattice)lattice).arrayIndex(coord));
        }
    };
}
