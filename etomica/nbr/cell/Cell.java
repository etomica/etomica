package etomica.nbr.cell;

import etomica.api.IAtom;
import etomica.atom.AtomArrayList;
import etomica.lattice.AbstractLattice;
import etomica.lattice.RectangularLattice;
import etomica.lattice.SiteFactory;
import etomica.util.Debug;

/**
 * Site used to form array of cells for cell-based neighbor listing.  Each
 * cell is capable of holding lists of atoms that are in them.
 */
public class Cell implements java.io.Serializable {

    public Cell(int latticeArrayIndex) {
        this.latticeArrayIndex = latticeArrayIndex;
    }
    
    public AtomArrayList occupants() {
        if (occupants.getAtomCount() > 1 && occupants.getAtom(0) == occupants.getAtom(1)) {
            throw new RuntimeException("oops");
        }
        return occupants;
    }
    
    public void addAtom(IAtom atom) {
        if (Debug.ON && (occupants.getAtomCount() > 0 && occupants.getAtom(0) == atom)) {
            throw new RuntimeException("oops "+atom);
        }
        occupants.add(atom);
    }
    
    public void removeAtom(IAtom atom) {
        occupants.removeAndReplace(occupants.indexOf(atom));
    }
    
    public int getLatticeArrayIndex() {
        return latticeArrayIndex;
    }
    
    private final AtomArrayList occupants = new AtomArrayList(1);
    final int latticeArrayIndex;//identifies site in lattice

    private static final long serialVersionUID = 1L;
    public static final SiteFactory FACTORY = new CellFactory();
    
    public static class CellFactory implements SiteFactory, java.io.Serializable {
        private static final long serialVersionUID = 1L;

        public Object makeSite(AbstractLattice lattice, int[] coord) {
            return new Cell(((RectangularLattice)lattice).arrayIndex(coord));
        }
    };
}
