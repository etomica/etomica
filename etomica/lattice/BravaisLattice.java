package etomica.lattice;

import etomica.Space;

/**
 * Arbitrary-dimension Bravais Lattice, in which the sites are instances of 
 * Space.Vector, with positions given as linear combinations of a set of
 * primitive vectors.
 */

public class BravaisLattice implements SpaceLattice {

    public BravaisLattice(Primitive primitive) {
//        super(primitive.space.D(), new BravaisSiteFactory());
        this.primitive = primitive;
    }

    public int D() {
        return getSpace().D();
    }
    
    public Space getSpace() {
        return primitive.space;
    }
    /**
     * Calculates and returns a new vector that is the spatial position given
     * by adding together the primitive vectors, each multiplied by the corresponding
     * integer index given by the array argument.
     */
    public Object site(int[] index) {
        if(index.length != getSpace().D()) throw new IllegalArgumentException("index given to site method of lattice must have number of elements equal to dimension of lattice");
        Space.Vector vector = getSpace().makeVector();
        for(int i=0; i<index.length; i++) {
            vector.PEa1Tv1(index[i], primitive.latticeVectors[i]);
        }
        return vector;
    }

    

    /**
     * Sets the primitive for this lattice to the one given, and
     * updates the site positions.
     */
    public void setPrimitive(Primitive primitive) {
        this.primitive = primitive;
//        setSize(size);//rebuilds sites
    }
    
    /**
     * Returns the primitive object used to construct this lattice. Note that if
     * the primitive is modified, changes will not be reflected in this lattice
     * until the update() method is called.
     */
    public Primitive getPrimitive() {
        return primitive;
    }

//    //not carefully implemented
//    public Space.Vector nearestSite(Space.Vector r) {
//        Space.Vector[] coords = (Space.Vector[])sites();
//        double r2min = Double.MAX_VALUE;
//        Space.Vector site = null;
//        for(int i=coords.length-1; i>=0; i--) {
//            double r2 = coords[i].Mv1Squared(r);
//            if (r2 < r2min) {
//                r2min = r2;
//                site = coords[i];
//            }
//        }
//        return site;
//    }

//    /**
//     * Returns the event manager that registers listeners and notifies them of
//     * events indicating changes in this lattice. Part of AbstractLattice
//     * interface.
//     */
//    public SimulationEventManager eventManager() {
//        return eventManager;
//    }
//
//    /**
//     * Returns "BravaisLattice" plus the array size, e.g., BravaisLattice{20,20,10}.
//     */
//    public String toString() {
//        return "BravaisLattice"+Arrays.toString(size);
//    }
//
//    private static class BravaisSiteFactory implements SiteFactory {
//        public Object makeSite(AbstractLattice lattice, int[] index) {
//            return ((BravaisLattice)lattice).getPrimitive().position(index);
//        }
//    }

    private Primitive primitive;
//    public final SimulationEventManager eventManager = new SimulationEventManager();
//    private final LatticeEvent rebuildEvent = new LatticeEvent(this,
//            LatticeEvent.REBUILD);
//    private final LatticeEvent allSiteEvent = new LatticeEvent(this,
//            LatticeEvent.ALL_SITE);
//    private final LatticeEvent resetNbrEvent = new LatticeEvent(this,
//            LatticeEvent.RESET_NBRS);
    
}//end of BravaisLattice
