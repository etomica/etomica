package simulate;

public interface Lattice {
    
    public int D();                      //dimension (1D, 2D, etc) of the lattice
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
        public Iterator adjacentIterator();   //iterator for all sites adjacent to the site
        public Lattice lattice();             //lattice on which the site resides
        public boolean isAdjacentTo(Site s);  //test for adjacency of the site to another site
        public Coordinate coordinate();       //coordinate of the site
        
        public interface Habitable extends Site {
            public Occupant[] occupants();
        }
        
        public interface Iterator {
            public boolean hasNext();         //true if another iterate in list
            public void reset();              //put iterator back in beginning state
            public Site first();              //return first element of list without affecting state of iterator
            public Site next();               //advance iterator and return next element
            public void allSites(Action act); //method to perform an action sequentially to every site given by iterator
        }
        
        public static final class Linker {    //class for making linked lists of sites
            public final Lattice.Site site;
            public final Linker next;
            public Linker(Linker sl, Lattice.Site s) {next = sl; site = s;}
        }
        
        public interface Action {
            public void action(Site s);
        }
    }
}