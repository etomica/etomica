package etomica;

/**
 * Determines and keeps lists of neighbors uplist and downlist of a 
 * reference atom.  "up" and "down" distinguish neighbor by whether they
 * come before (down) or after (up) the reference atom when looping over
 * a given iterator.  Definition of "neighbor" is given by neighbor criterion
 * passed to setNeighbors method.
 */
public class NeighborManager {
    
    private final AtomList neighborList;
    private final Atom atom;
 //   public final AtomLinker.Tab tab;
    public final AtomLinker tab;
    
    public NeighborManager(Atom s) {
        atom = s;
        neighborList = new AtomList();
        tab = new AtomLinker.Tab();
 //       tab = new AtomLinker(s);
        neighborList.add(tab);
    }
    public NeighborManager(Atom s, AtomList list, Criterion criterion) {
        this(s);
        setupNeighbors(list, criterion);
    }

    public Atom atom() {return atom;}
    public int neighborCount() {return neighborList.size();}
    public boolean isNeighbor(Atom s) {
        return (neighborList.contains(s));
    }
    
    public AtomList neighbors() {return neighborList;}
    
    public void clearAll() {
        neighborList.clear();
        neighborList.add(tab);
    }

    public void setupNeighbors(AtomList list, Criterion criterion) {  //set up neighbors according to given criterion
        clearAll();
        AtomIteratorList iterator = new AtomIteratorList(list);
        boolean down = true;
        while(iterator.hasNext()) {              //begin outer loop
            Atom s = iterator.next();
            if(s == atom) {down = false;}     //subsequent neighbors go in up-list
            else if(criterion.areNeighbors(s,atom)) {
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
            public boolean areNeighbors(Atom s1, Atom s2);
             
            /**
             * Criterion that defines all sites on the lattice to be neighbors of each other.
             */
            public class All implements Criterion {
                public boolean areNeighbors(Atom s1, Atom s2) {return s1 != s2;}
            }//end of Criterion.All
        }//end of NeighborManager.Criterion
                    
}//end of NeighborManager
            
