package etomica.lattice;

public class Site implements java.io.Serializable {
    
    private AbstractLattice lattice;
    private SiteIterator.Neighbor neighborIterator;
    private AbstractLattice.Coordinate coordinate;
    
    /**
     * Creates a site having the given parent lattice and coordinate.
     * Neighbor-iterator is set to an empty SiteIterator.Neighbor class, and may be subsequently filled or reassigned.
     */
    public Site(AbstractLattice parent, AbstractLattice.Coordinate coord) {
        this(parent, new SiteIterator.Neighbor(), coord);
    }
    public Site(AbstractLattice parent, SiteIterator.Neighbor iterator, AbstractLattice.Coordinate coord) {
        lattice = parent;
        neighborIterator = iterator;
        this.coordinate = coord;
    }
    
    public AbstractLattice lattice() {return lattice;}  //returns the (top-level) lattice on which this site resides
    public SiteIterator.Neighbor neighborIterator() {return neighborIterator;}
    public void setNeighborIterator(SiteIterator.Neighbor iterator) {neighborIterator = iterator;}
    public AbstractLattice.Coordinate coordinate() {return coordinate;}
    public void setCoordinate(AbstractLattice.Coordinate coord) {coordinate = coord;}
    /**
     * Test for adjacency of the site to another site
     */
    public boolean isNeighborOf(Site s) {
        neighborIterator.reset();
        while(neighborIterator.hasNext()) {
            if(s == neighborIterator.next()) {return true;}
        }
        return false;
    }
    
    public String toString() {return coordinate.toString();}
    
    ///////// end of Site methods and fields ///////////
    
    public static class Factory implements SiteFactory {
        public Site makeSite(AbstractLattice parent, SiteIterator.Neighbor iterator, AbstractLattice.Coordinate coord) {
            return new Site(parent, iterator, coord);
        }
    }
        
//    public interface Habitable extends Site {
//        public LatticeAbstract.Occupant[] occupants();
//    }
        
    public static final class Linker {    //class for making linked lists of sites
        public final Site site;
        public Linker next;
        public Linker(Site s, Linker sl) {next = sl; site = s;}
    }
}