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
    
    private final AtomList neighborList;
    private final Site site;
 //   public final AtomLinker.Tab tab;
    public final AtomLinker tab;
    
    public NeighborManager(Site s) {
        site = s;
        neighborList = new AtomList();
 //       tab = new AtomLinker.Tab();
        tab = new AtomLinker(s);
        neighborList.add(tab);
    }
    public NeighborManager(Site s, AtomList list, Criterion criterion) {
        this(s);
        setupNeighbors(list, criterion);
    }

    public Site site() {return site;}
    public int neighborCount() {return neighborList.size();}
    public boolean isNeighbor(Site s) {
        return (neighborList.contains(s));
    }
    
    /**
     * Returns the first neighbor that would be encountered when proceeding
     * up the iterator used to construct the neighbor lists.
     */
    public Atom firstNeighbor() {return neighborList.getFirst();}
    /**
     * Returns the last neighbor that would be encountered when proceeding
     * up the iterator used to construct the neighbor lists.
     */
    public Atom lastNeighbor() {return neighborList.getLast();}
    
    public AtomList neighbors() {return neighborList;}
    
    public void clearAll() {
        neighborList.clear();
    }

    public void setupNeighbors(AtomList list, Criterion criterion) {  //set up neighbors according to given criterion
        AtomIteratorList iterator = new AtomIteratorList(list);
        iterator.reset();
        boolean down = true;
        while(iterator.hasNext()) {              //begin outer loop
            Site s = (Site)iterator.next();
            if(s == site) {down = false;}     //subsequent neighbors go in up-list
            else if(criterion.areNeighbors(s,site)) {
                if(down) neighborList.addBefore(s, tab);
                else     neighborList.addLast(s);
            }
        }
    }//end of setupNeighbors
            
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
        }//end of NeighborManager.Criterion
                    
}//end of NeighborManager
            
