/*
 * History
 * Created on Nov 23, 2004 by kofke
 */
package etomica.nbr.site;

import etomica.atom.Atom;
import etomica.lattice.AbstractLattice;
import etomica.lattice.RectangularLattice;
import etomica.lattice.SiteFactory;

/**
 * Site used to form array of cells for cell-based neighbor listing.  Each
 * cell is capable of holding lists of atoms that are in them.
 */


public class AtomSite {

    public AtomSite(int latticeArrayIndex) {
        this.latticeArrayIndex = latticeArrayIndex;
    }
    
    public Atom getAtom() {return atom;}
    
    public void setAtom(Atom atom) {
        this.atom = atom;
        if (atom == null) return;
        AtomSequencerSite seq = (AtomSequencerSite)atom.seq;
        if(this == seq.site) return;
        if(seq.site != null) seq.site.setAtom(null);
        seq.site= this;
    }
    
    public int getLatticeArrayIndex() {
        return latticeArrayIndex;
    }
    
    private Atom atom;
    final int latticeArrayIndex;//identifies site in lattice

    public static final SiteFactory FACTORY = new SiteFactory() {
        public Object makeSite(AbstractLattice lattice, int[] coord) {
            return new AtomSite(((RectangularLattice)lattice).arrayIndex(coord));
        }
    };
}
