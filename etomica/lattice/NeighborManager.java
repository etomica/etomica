package etomica.lattice;
import etomica.*;

/**
 * Determines and keeps lists of neighbors uplist and downlist of a 
 * reference site.  "up" and "down" distinguish neighbor by whether they
 * come before (down) or after (up) the reference site when looping over
 * a given iterator.  Definition of "neighbor" is given by neighbor criterion
 * passed to setNeighbors method.
 */
public class NeighborManager {
    
    private final AtomList upList;
    private final AtomList dnList;
    private final Site site;
    
    public NeighborManager(Site s) {
        site = s;
        upList = new AtomList();
        dnList = new AtomList();
    }
    public NeighborManager(Site s, AtomIterator iterator, Criterion criterion) {
        this(s);
        setupNeighbors(iterator, criterion);
    }

    public Site site() {return site;}
    public int neighborCount() {return upList.size() + dnList.size();}
    public boolean isNeighbor(Site s) {
        return (upList.contains(s) || dnList.contains(s));
    }
    
    /**
     * Returns the first neighbor that would be encountered when proceeding
     * up the iterator used to construct the neighbor lists.
     */
    public Site firstNeighbor() {return (Site)dnList.getLast();}
    /**
     * Returns the last neighbor that would be encountered when proceeding
     * up the iterator used to construct the neighbor lists.
     */
    public Site lastNeighbor() {return (Site)upList.getLast();}
    
    public AtomList upNeighbors() {return upList;}
    public AtomList downNeighbors() {return dnList;}
    
    public void clearAll() {
        upList.clear();
        dnList.clear();
    }

    public void setupNeighbors(AtomIterator iterator, Criterion criterion) {  //set up neighbors according to given criterion
        iterator.reset();
        boolean down = true;
        while(iterator.hasNext()) {              //begin outer loop
            Site s = (Site)iterator.next();
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
            
