package etomica.lattice;

import etomica.*;

public class Site extends Atom {
    
    private AbstractLattice.Coordinate coordinate;
    private final NeighborManager neighborManager = new NeighborManager(this);
    public int index;//temporary -- put in sequencer
    
    /**
     * Creates a site having the given parent lattice and coordinate.
     */
    public Site(Space space, AtomType type, AtomTreeNode.Factory nodeFactory) {
        super(space, type, nodeFactory, IteratorFactorySimple.INSTANCE.atomSequencerFactory());
    }
        
    public AbstractLattice.Coordinate coordinate() {return coordinate;}
    public void setCoordinate(AbstractLattice.Coordinate coord) {coordinate = coord;}

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