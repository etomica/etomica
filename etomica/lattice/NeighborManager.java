package etomica.lattice;
import etomica.*;
import etomica.AtomLinker.Tab;

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
    public /*final*/ AtomLinker tab;
    public final static int NEIGHBOR_TAB = Tab.requestTabType(); 

    public NeighborManager(Site s) {
        site = s;
        neighborList = new AtomList();
 //       tab = selfNeighbor ? new AtomLinker(s) : AtomLinker.newTab();
 //       neighborList.add(tab);
    }
        
    public Site site() {return site;}
    public int neighborCount() {return neighborList.size();}
    public boolean isNeighbor(Site s) {
        return (neighborList.contains(s));
    }
    
    public AtomList neighbors() {return neighborList;}
    
    public void clearAll() {
        neighborList.clear();
        if(tab != null) neighborList.add(tab);
    }

    public void setupNeighbors(AtomList list, Criterion criterion) {  //set up neighbors according to given criterion
  //      clearAll();
        neighborList.clear();
        AtomIteratorList iterator = new AtomIteratorList(list);
      if(criterion.areNeighbors(site, site)) {//make self neighbor
          if(tab == null || tab instanceof AtomLinker.Tab) tab = new AtomLinker(site);
      } else {//not self neighbor
          if(tab == null || !(tab instanceof AtomLinker.Tab)) tab = AtomLinker.newTab(neighborList, NEIGHBOR_TAB);
      }
      neighborList.add(tab);
      boolean down = true;
        while(iterator.hasNext()) {              //begin outer loop
            Site s = (Site)iterator.nextAtom();
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
            
