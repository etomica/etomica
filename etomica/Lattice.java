package simulate;

public interface Lattice {
        
    public int siteCount();              //total number of sites on lattice
    public int coordinationNumber();     //number of neighbors of each site
    public Site site(Coordinate coord);  //obtain a site given some specification in coordinate
    public Site randomSite();            //get a random site
    public Site.Iterator iterator();     //iterator for all sites in lattice
    
    public interface Coordinate {}
    
    public interface Occupant {
        public Site site();
    }
    
    public interface Site {
        public Iterator neighborIterator();
        public Lattice parentLattice();
        
        public interface Iterator {
            public boolean hasNext();
            public void reset();
            public Site first();
            public Site next();
        }
    }
}