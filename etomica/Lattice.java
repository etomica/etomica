package simulate;

public interface Lattice {
    
    public static boolean PERIODIC = true;
    public static boolean NOT_PERIODIC = false;
    
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
        public NeighborIterator adjacentIterator();   //iterator for all sites adjacent to the site
        public Lattice lattice();             //lattice on which the site resides
        public boolean isAdjacentTo(Site s);  //test for adjacency of the site to another site
        public Coordinate coordinate();       //coordinate of the site
        
        public interface Habitable extends Site {
            public Occupant[] occupants();
        }
        
        public interface Iterator extends Cloneable {
            public boolean hasNext();         //true if another iterate in list
            public void reset();              //put iterator back in beginning state
            public Site first();              //return first element of list without affecting state of iterator
            public Site next();               //advance iterator and return next element
            public void allSites(Action act); //method to perform an action sequentially to every site given by iterator
            public Object clone();            //intended to create a shallow copy of the iterator, which permits iteration over the same set of sites independently of (and without interfering with) the original
        }

        //Iterator to generate adjacent sites on LatticeBravais
        public static class NeighborIterator implements Iterator {
            private Lattice.Site.Linker firstUp, firstDown, nextLink;
            private boolean doBoth = false;
            private int count;
            private final Site site;
            public NeighborIterator(Site s) {site = s;}
            public Site site() {return site;}
            public int neighborCount() {return count;}
            public boolean hasNext() {return nextLink != null;}
            public void reset() {nextLink = firstUp; doBoth = true;}
            public void resetUp() {nextLink = firstUp; doBoth = false;} //put iterator in state ready to traverse up-neighbors
            public void resetDown() {nextLink = firstDown; doBoth = false;} //put iterator in state ready to traverse down-neighbors
            public Lattice.Site first() {return firstUp.site;}
            public Lattice.Site firstUp() {return firstUp.site;} //return first element of up-list without affecting state of iterator
            public Lattice.Site firstDown() {return firstDown.site;} //return first element of down-list without affecting state of iterator
            public Lattice.Site next() {
                Lattice.Site nextSite = nextLink.site;
                nextLink = nextLink.next;
                if(nextLink==null && doBoth==true) {nextLink = firstDown; doBoth = false;}
                return nextSite;
            }
            public void allSitesUp(Lattice.Site.Action act) {  //method to perform an action sequentially to every up-site given by iterator
                for(Lattice.Site.Linker l=firstUp; l!=null; l=l.next) {act.action(l.site);}
            }
            public void allSitesDown(Lattice.Site.Action act) { //method to perform an action sequentially to every down-site given by iterator
                for(Lattice.Site.Linker l=firstDown; l!=null; l=l.next) {act.action(l.site);}
            }
            public void allSites(Lattice.Site.Action act) {
                for(Lattice.Site.Linker l=firstUp; l!=null; l=l.next) {act.action(l.site);}
                for(Lattice.Site.Linker l=firstDown; l!=null; l=l.next) {act.action(l.site);}
            }
            public void addUp(Lattice.Site s) {firstUp = new Lattice.Site.Linker(firstUp, s); count++;}
            public void addDown(Lattice.Site s) {firstDown = new Lattice.Site.Linker(firstDown, s); count++;}
            public void setNeighbors(Iterator iterator, Criterion criterion) {  //set up neighbors according to given criterion
                iterator.reset();
                boolean down = true;
                while(iterator.hasNext()) {              //begin outer loop
                    Site s = iterator.next();
                    if(s == site) {down = false;}     //subsequent neighbors go in up-list
                    else if(criterion.isNeighbor(s,site)) {
                        if(down) {addDown(s);}
                        else     {addUp(s);}
                    }
                }
            }
            public Object clone() { 
                try {return super.clone();} catch(CloneNotSupportedException e) {return null;}
            }
                
            
            public interface Criterion {
                public boolean isNeighbor(Site s1, Site s2);
                
                public class All implements Criterion {
                    public boolean isNeighbor(Site s1, Site s2) {return s1 != s2;}
                }
            }
        }
                
        public static final class Linker {    //class for making linked lists of sites
            public final Lattice.Site site;
            public Linker next;
            public Linker(Linker sl, Lattice.Site s) {next = sl; site = s;}
        }
        
        public interface Action {
            public void action(Site s);
        }
    }
}