/*
 * History
 * Created on Nov 23, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomList;
import etomica.lattice.AbstractLattice;
import etomica.lattice.RectangularLattice;
import etomica.lattice.SiteFactory;

/**
 * Site used to form array of cells for cell-based neighbor listing.  Each
 * cell is capable of holding lists of atoms that are in them.
 */


public class Cell implements java.io.Serializable {

    public Cell(int latticeArrayIndex) {
        this.latticeArrayIndex = latticeArrayIndex;
    }
    
    public AtomArrayList occupants() {return occupants;}
    
    public void addAtom(Atom atom) {
        occupants.add(atom);
    }
    
    public int getLatticeArrayIndex() {
        return latticeArrayIndex;
    }
    
    private final AtomArrayList occupants = new AtomArrayList();
    final int latticeArrayIndex;//identifies site in lattice

    public static final SiteFactory FACTORY = new CellFactory();
    
    public static class CellFactory implements SiteFactory, java.io.Serializable {
        public Object makeSite(AbstractLattice lattice, int[] coord) {
            return new Cell(((RectangularLattice)lattice).arrayIndex(coord));
        }
    };
}
