package etomica.lattice;

import etomica.SimulationEventManager;
import etomica.Space;
import etomica.utility.Arrays;

/**
 * Arbitrary-dimension Bravais Lattice, in which the sites are instances of 
 * Space.Vector, with positions given as linear combinations of a set of
 * primitive vectors.
 */

public class BravaisLattice extends SimpleLattice {

    /**
     * Constructor should be invoked only within a build() method of a Factory
     * class. The build method handles the construction of the tree structure
     * under this instance, which forms the lattice.
     */
    public BravaisLattice(Primitive primitive) {
        super(primitive.space.D(), new BravaisSiteFactory());
        this.primitive = primitive;
    }

    /**
     * Returns a site roughly at the center of the lattice.
     */
    public Object centerSite() {
        int[] halfDim = new int[D()];
        for (int i = 0; i < D(); i++)
            halfDim[i] = (int) size[i] / 2;
        return site(halfDim);
    }

    /**
     * Sets the size of the lattice (number of atoms in each direction) so that
     * it has a number of sites equal or greater than the given value n. Sets
     * lattice to be square (same number in all directions), and finds smallest
     * size that gives number of sites equal or exceeding the given value.
     */
    public void setSize(int n) {
        if (n < 0 || n >= Integer.MAX_VALUE)
            throw new IllegalArgumentException(
                    "Inappropriate size specified for lattice: " + n);
        int i = (int)Math.pow(n,1/D());
        while(i <= n) {
            if ((int)Math.pow(i,D()) >= n) break;
            i++;
        }
        setSize(primitive.space.makeArrayD(i));
    }

    /**
     * Performs same actions as setSize(int), then size of primitive vectors are
     * adjusted such that lattice will fit in the given dimensions. Assumes all
     * dimension values are equal (future development will address case of non-
     * square dimensions).
     */
    public void setSize(int n, double[] dimensions) {
        if (n < 0 || n >= Integer.MAX_VALUE)
            throw new IllegalArgumentException(
                    "Inappropriate size specified for lattice: " + n);
        int i = (int)Math.pow(n,1/D());
        while(i <= n) {
            if ((int)Math.pow(i,D()) >= n) break;
            i++;
        }
        //size primitive vectors to given dimensions
        double[] newLength = new double[D()];
        for (int j = 0; j < D(); j++)
            newLength[j] = dimensions[j] / (double) i;
        primitive.setSize(newLength);
        setSize(primitive.space.makeArrayD(i));
    }

    /**
     * Translates the lattice so that the first site is positioned at the
     * given point.
     */
    public void shiftFirstTo(Space.Vector r) {
        Space.Vector[] coords = (Space.Vector[])sites();
        for(int i=coords.length-1; i>=0; i--) {
            coords[i].ME(r);
        }
    }

    /**
     * Sets the primitive for this lattice to the one given, and
     * updates the site positions.
     */
    public void setPrimitive(Primitive primitive) {
        this.primitive = primitive;
        setSize(size);//rebuilds sites
    }
    
    /**
     * Returns the primitive object used to construct this lattice. Note that if
     * the primitive is modified, changes will not be reflected in this lattice
     * until the update() method is called.
     */
    public Primitive getPrimitive() {
        return primitive;
    }

    //not carefully implemented
    public Space.Vector nearestSite(Space.Vector r) {
        Space.Vector[] coords = (Space.Vector[])sites();
        double r2min = Double.MAX_VALUE;
        Space.Vector site = null;
        for(int i=coords.length-1; i>=0; i--) {
            double r2 = coords[i].Mv1Squared(r);
            if (r2 < r2min) {
                r2min = r2;
                site = coords[i];
            }
        }
        return site;
    }

    /**
     * Returns the event manager that registers listeners and notifies them of
     * events indicating changes in this lattice. Part of AbstractLattice
     * interface.
     */
    public SimulationEventManager eventManager() {
        return eventManager;
    }

    /**
     * Returns "BravaisLattice" plus the array size, e.g., BravaisLattice{20,20,10}.
     */
    public String toString() {
        return "BravaisLattice"+Arrays.toString(size);
    }

    private static class BravaisSiteFactory implements SiteFactory {
        public Object makeSite(AbstractLattice lattice, int[] index) {
            return ((BravaisLattice)lattice).getPrimitive().position(index);
        }
    }

    private Primitive primitive;
    public final SimulationEventManager eventManager = new SimulationEventManager();
    private final LatticeEvent rebuildEvent = new LatticeEvent(this,
            LatticeEvent.REBUILD);
    private final LatticeEvent allSiteEvent = new LatticeEvent(this,
            LatticeEvent.ALL_SITE);
    private final LatticeEvent resetNbrEvent = new LatticeEvent(this,
            LatticeEvent.RESET_NBRS);
    
}//end of BravaisLattice
