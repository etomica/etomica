package etomica.lattice;

import etomica.*;

public class Site extends Atom {
    
    private final NeighborManager neighborManager = new NeighborManager(this);
    
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
        return node.index() < anotherSite.node.index();
    }
    
    public int D() {return coord.position().D();}
    
    
    public int[] latticeCoordinate() {
        int D = D();
        int[] idx = new int[D];
        AtomTreeNode node = this.node;
        for(int i=D-1; i>=0; i--, node=node.parentNode()) {
            idx[i] = node.index();
        }
        return idx;
    }//end of latticeCoordinate
        
public static class Factory extends AtomFactory {
    
    AtomType atomType;
    
    public Factory(Space space) {
        super(space);
        setType(new AtomType(this));//default
    }
    
    public void setType(AtomType t) {atomType = t;}
    public AtomType type() {return atomType;}
    
    /**
     * Builds a single site.
     */
    protected Atom build() {
        return new Site(space, atomType);
    }
    
}//end of AtomFactorySite
        
}//end of Site