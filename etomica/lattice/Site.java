package etomica.lattice;

import etomica.*;

public class Site extends Atom {
    
    private final NeighborManager neighborManager = new NeighborManager(this);
    public int index;//temporary -- put in sequencer
    
    /**
     * Creates a site having the given parent lattice and coordinate.
     */
    public Site(Space space, AtomType type) {
        super(space, type, AtomTreeNodeLeaf.FACTORY, IteratorFactorySimple.INSTANCE.atomSequencerFactory());
    }
        
    public NeighborManager neighborManager() {return neighborManager;}
    /**
     * Test for adjacency of the site to another site
     */
    public boolean isNeighborOf(Site s) {
        return neighborManager.isNeighbor(s);
    }
    
    //temporary -- put in sequencer
    public boolean preceeds(Site anotherSite) {
        return index < anotherSite.index;
    }
        
}//end of Site