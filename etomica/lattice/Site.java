package etomica.lattice;

public class Site implements java.io.Serializable {
    
    private final AbstractLattice lattice;
    private AbstractLattice.Coordinate coordinate;
    private final NeighborManager neighborManager = new NeighborManager(this);
    
    /**
     * Creates a site having the given parent lattice and coordinate.
     */
    public Site(AbstractLattice parent, AbstractLattice.Coordinate coord) {
        lattice = parent;
        this.coordinate = coord;
    }
    
    public AbstractLattice lattice() {return lattice;}  //returns the (top-level) lattice on which this site resides
    public AbstractLattice.Coordinate coordinate() {return coordinate;}
    public void setCoordinate(AbstractLattice.Coordinate coord) {coordinate = coord;}

    public NeighborManager neighborManager() {return neighborManager;}
    /**
     * Test for adjacency of the site to another site
     */
    public boolean isNeighborOf(Site s) {
        return neighborManager.isNeighbor(s);
    }
    
    public String toString() {return coordinate.toString();}
    
    ///////// end of Site methods and fields ///////////
    
    public static class Factory implements SiteFactory {
        public Site makeSite(AbstractLattice parent, AbstractLattice.Coordinate coord) {
            return new Site(parent, coord);
        }
    }
        
}//end of Site