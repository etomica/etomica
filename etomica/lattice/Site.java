package etomica.lattice;

import etomica.*;

public class Site extends Atom {
    
    private final NeighborManager neighborManager = new NeighborManager(this);
    
    /**
     * Creates a site having the given parent lattice and coordinate.
     * Sets node such that site is a leaf, not a group.
     */
    public Site(Space space, AtomType type, AtomTreeNodeGroup parent) {
        this(space, type, AtomTreeNodeLeaf.FACTORY, parent);
    }
    /**
     * Creates a site having the given parent lattice and coordinate.
     */
    public Site(Space space, AtomType type, AtomTreeNode.Factory nodeFactory, AtomTreeNodeGroup parent) {
        super(space, type, nodeFactory, 
            IteratorFactorySimple.INSTANCE.simpleSequencerFactory(), parent);
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
    
    /**
     * Returns the spatial dimension of the lattice on
     * which this site resides.
     */
    public int D() {return coord.position().D();}
    
    /**
     * Returns the lattice coordinate of this site, defined such
     * that upon passing the coordinate to the site(int[]) method
     * of the parent lattice, this site would be returned.
     */
    public int[] latticeCoordinate() {
        int D = D();
        int[] idx = new int[D];
        AtomTreeNode node = this.node;
        for(int i=D-1; i>=0; i--, node=node.parentNode()) {
            idx[i] = node.index();
        }
        return idx;
    }//end of latticeCoordinate

///////////////////////////////////////////////////////////////////////////

/**
 * Factory that constructs a simple Site.
 */
public static class Factory extends AtomFactory {
    
    AtomType atomType;
    
    public Factory(Simulation sim) {
        super(sim);
        setType(new AtomType(this));//default
    }
    
    public boolean isGroupFactory() {return false;}
    
    public void setType(AtomType t) {atomType = t;}
    public AtomType type() {return atomType;}
    
    /**
     * Builds a single site.
     */
    protected Atom build(AtomTreeNodeGroup parent) {
        return new Site(parentSimulation().space, atomType, parent);
    }

    public Atom build(Atom atom) {
        if(atom.creator() != this) throw new IllegalArgumentException("Site.factory.build(Atom) attempted using an atom that was not built by the same factory");
        return atom;
    }
    
}//end of AtomFactorySite
        
}//end of Site