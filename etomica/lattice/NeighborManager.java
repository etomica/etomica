package etomica.lattice;

/**
 * Determines and keeps lists of neighbors uplist and downlist of a 
 * reference site.  "up" and "down" distinguish neighbor by whether they
 * come before (down) or after (up) the reference site when looping over
 * a given iterator.
 */
public class NeighborManager {
    
    private final SiteList upList;
    private final SiteList dnList;
    private final Site site;
    
    public NeighborManager(Site s) {
        site = s;
        upList = new SiteList();
        dnList = new SiteList();
    }
    public NeighborManager(Site s, SiteIterator iterator, Criterion criterion) {
        this(s);
        setNeighbors(iterator, criterion);
    }

    public Site site() {return site;}
    public int neighborCount() {return upList.size() + dnList.size();}
    
    public SiteList upNeighbors() {return upList;}
    public SiteList downNeighbors() {return dnList;}
    
    public void clearAll() {
        upList.clear();
        dnList.clear();
    }

    public void setNeighbors(SiteIterator iterator, Criterion criterion) {  //set up neighbors according to given criterion
        iterator.reset();
        boolean down = true;
        while(iterator.hasNext()) {              //begin outer loop
            Site s = iterator.next();
            if(s == site) {down = false;}     //subsequent neighbors go in up-list
            else if(criterion.areNeighbors(s,site)) {
                if(down) {dnList.addFirst(s);}
                else     {upList.addLast(s);}
            }
        }
    }//end of SiteIterator.Neighbor.setNeighbors
            
        /**
         * Defines a criterion for specifying whether two sites on a lattice are to be designated as neighbors of each other.
         * Used by the NeighborManager to construct lists of all "neighbors" of a given site.
         * The "neighbors" storted in the lists are those for which the areNeighbors method of this class returns 
         * <code>true</code>.
         */
        public interface Criterion {
            public boolean areNeighbors(Site s1, Site s2);
             
            /**
             * Criterion that defines all sites on the lattice to be neighbors of each other.
             */
            public class All implements Criterion {
                public boolean areNeighbors(Site s1, Site s2) {return s1 != s2;}
            }//end of Criterion.All
        }//end of Neighbor.Criterion
                    
}//end of NeighborManager
            
