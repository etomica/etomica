package simulate.lattice;
import simulate.Space;

public interface AbstractLattice extends java.io.Serializable {
    
    public static boolean PERIODIC = true;
    public static boolean NOT_PERIODIC = false;
    
    public int D();                      //dimension (1D, 2D, etc) of the lattice
    public int siteCount();              //total number of sites on lattice
//    public int coordinationNumber();     //number of neighbors of each site
    public Site site(Coordinate coord);  //obtain a site given some specification in coordinate
    public Site randomSite();            //get a random site
    public SiteIterator iterator();     //iterator for all sites in lattice
    
    public interface Coordinate {}
    public interface PositionCoordinate extends Coordinate {
        public Space.Vector position();
    }
    
    public interface Occupant {
        public Site site();
    }
    
}