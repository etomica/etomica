package etomica.lattice;
import etomica.*;

public class LatticeCubicFcc extends Atom {
    
   private double latticeConstant;
    
   private AtomList siteList;
   private int[] dimensions;
   private final int D = 3;
   
   /**
    * Constructor should be invoked only within a build() method of the
    * Factory class.  The build method handles the construction of the
    * substructure under this instance, which forms the lattice.
    */
   private LatticeCubicFcc(Space space, AtomType type, int[] dimensions, AtomTreeNodeGroup parent) {
        super(space, type, parent);
        if(space.D() != 3) throw IllegalArgumentException("Error in LatticeCubicFcc constructor:  Given space is not 3D");
        idx = new int[D];
        this.dimensions = new int[dimensions.length];
        System.arraycopy(dimensions, 0, this.dimensions, 0, dimensions.length);
   }
   
    /**
     * Creates an fcc lattice with the given number of sites in a box of the given size,
     * using a simple site factory (Site.Factory).
     */
    public static LatticeCubicFcc makeLattice(Simulation sim, int nSites, double boxSize) {
        return makeLattice(sim, getSize(nSites), getLatticeConstant(nSites,boxSize), new Site.Factory()); 
    }
    /**
     * Creates an fcc lattice with the given number of sites in a box of the given size,
     * using the given site factory.
     */
    public static LatticeCubicFcc makeLattice(Simulation sim, AtomFactory siteFactory, 
                                                int nSites, double boxSize) {
        return makeLattice(sim, getSize(nSites), getLatticeConstant(nSites,boxSize), factory); 
    }
    /**
     * Creates an fcc lattice with 4*dimensions[1]*dimensions[2]*dimensions[3] sites, using the given lattice constant
     * and a simple site factory.
     */
    public static LatticeCubicFcc makeLattice(Simulation sim, int[] dimensions, double latticeConstant) {
        return makeLattice(sim, new Site.Factory(), dimensions, latticeConstant);
    }
    
    public static LatticeCubicFcc makeLattice(Simulation sim, AtomFactory siteFactory, 
                                                int[] dimensions, double latticeConstant) {
        return (LatticeCubicFcc)new Factory(sim, siteFactory, dimensions, latticeConstant).makeAtom();
    }
    
    /**
     * Determines the size of the cubic lattice consistent with an fcc lattice having
     * the given number of sites.  
     * Returns the smallest k such that 4*k^3 is greater than or equal to n.
     */
    private static int[] getSize(int n) {
        if(n <= 0) return new int[] {0, 0, 0};
        int nLat = 0;
        int k = 0;
        while(n > nLat) {
            k++;
            nLat = 4*k*k*k;
        }
        return new int[] {k, k, k};
    }
    
    private static double getLatticeConstant(int nSites, double boxSize) {
        int n = getSize(nSites)[0];
        return boxSize/(double)n;
    }

    
    public int[] dimensions() {return dimensions;}
    
    public void setLatticeConstant(double a) {
        latticeConstant = a;
        update();
    }
    
    
    /**
     * Translates the lattice so that the first site is positioned at the origin,
     */
    public void shiftFirstToOrigin() {
        Space.Vector shift = Space.makeVector(D);
        shift.ME(siteList.getFirst().coord.position());
        this.coord.translateBy(shift);
    }
       
    //not carefully implemented
    public Atom nearestSite(Space.Vector r) {
        throw new RuntimeException("method LatticeCubicFcc.nearestSite not yet implemented");
    }
    
    /**
     * Returns the site corresponding to the given index.
     */
     //needs work to improve probable inefficiency with list.get
     //and to handle index values out of size of lattice
    public Site site(int[] idx) {
        Atom site = this;
        int i=0;
        do site = ((AtomTreeNodeGroup)site.node).childList.get(idx[i++]);
        while(!site.node.isLeaf() && i<idx.length);
        return (Site)site;
    }
    
    /**
     * Returns a list of all the Bravais lattice sites for this lattice.
     * Sites making up the lattice basis are not returned individually
     * if the basis is multiatomic.
     */
    public AtomList siteList() {return siteList;}
    
    /**
     * Returns the spatial dimension of this lattice.
     */
    public int D() {return D;}
    
    /**
     * Overrides the superclass method to return a simple string.
     * Terminates the chain of parentGroup calls that define the signature
     * of an atom in this lattice.
     */
    public String signature() {return "LatticeCubicFcc";}
    
    /**
     * Causes all coordinates to update their position vectors, and
     * notifies any observers that a change has occurred.
     */
    public void update() {
        type.creator().getConfiguration().initializePositions(this);
        notifyObservers();
    }//end of update
    
    public void setupNeighbors(NeighborManager.Criterion criterion) {
        AtomIteratorList iterator = new AtomIteratorList(siteList);
        iterator.reset();
        while(iterator.hasNext()) {
            Site site = (Site)iterator.next();
            site.neighborManager().setupNeighbors(siteList, criterion);
        }
    }//end of setupNeighbors
    
    public Space.Vector[] positions() {
        int n = siteCount();
        Space.Vector[] r = new Space.Vector[n];
        for(int i=0; i<n; i++) {r[i] = new Space3D.Vector();}
        SiteIterator iteratorsites = iterator();
        iteratorsites.reset();
        int i = 0;
        while (iteratorsites.hasNext()){
            Site site = iteratorsites.next();
            r[i++].E(((AbstractLattice.PositionCoordinate)site.coordinate()).position());
        }
        return r;
    }//end of positions


//////////////////////////////////////////////////////////////////////////////////

/**
 * Criterion that identifies neighbors of a site on a Bravais lattice as
 * those sites directly adjacent to it, where "adjacent" means that exactly
 * one index differs by +/- 1.  For example, in 2D the neighbors
 * of site (1,2) would be {(0,2),(1,1),(1,3),(2,2)}.  If periodic is 
 * <code>true</code>, then neighbors are indicated also if sites are adjacent
 * upon application of toroidal peroidic boundaries (thus, for example, site 
 * (0,1) would be a neighbor of site (3,1) if the lattice dimensions are (4,n)).
 * Default condition is such that periodic boundaries are applied.
 */

public static class AdjacencyCriterion implements NeighborManager.Criterion {
    
    private boolean periodic = true;
    private BravaisLattice lattice;
    
    public AdjacencyCriterion(BravaisLattice lattice) {
        this.lattice = lattice;
    }
    
    /**
     * Indicates if periodicity is applied in defining adjacency.
     */
    public boolean isPeriodic() {return periodic;}
    /**
     * Sets whether periodicity is applied in defining adjacency.
     * Default value is true.
     */
    public void setPeriodic(boolean b) {periodic = b;}
    
    /**
     * Returns true if the given sites are adjacent to each other,
     * considering periodicity as set previously.
     */
    public boolean areNeighbors(Site s1, Site s2) {
        int[] dim = lattice.dimensions();
        int[] idx1 = s1.latticeCoordinate();
        int[] idx2 = s2.latticeCoordinate();

        boolean sameTillNow = true;
        for(int i=0; i<dim.length; i++) {
            int imax = dim[i] - 1;
            int diff = idx1[i] - idx2[i];
            if(periodic && (diff == imax || diff == -imax) ) diff = +1;
            switch(diff) {
                case 0: break;
                case +1:
                case -1: if(sameTillNow) sameTillNow = false;//found first difference
                            else return false;//more than one differs
                            break;
                default: return false;//found one that differs by more than +/- 1
            }//end switch
        }//end for
        return !sameTillNow;//returns false if all indexes were the same, true if one differed by 1
    }//end areNeighbors
}//end AdjacencyCriterion

//////////////////////////////////////////////////////////////////////////////////
   
   
public static class Factory extends AtomFactory {
        
    private LatticeFactoryCubic subFactory;
    
    /**
     * Creates an fcc lattice with 4*dimensions[1]*dimensions[2]*dimensions[3] sites, 
     * using the given lattice constant and site factory.
     */
    public Factory(Simulation sim, AtomFactory siteFactory, int[] dimensions, double latticeConstant) {
        super(sim);
        this.latticeConstant = latticeConstant;
        subFactory = new LatticeFactoryCubic(sim, siteFactory, dimensions, latticeConstant);
        configuration = new Configuration4();
    }
    
    public Atom build(AtomTreeNodeGroup parent) {
        AtomGroup group = new LatticeCubicFcc(parentSimulation.space, groupType, parent);
        return build(group);
    }
    
    public Atom build(Atom atom) {
        if(!(atom instanceof LatticeCubicFcc)) throw new IllegalArgumentException("Error in LatticeFactoryCubicFcc.build(Atom):  Attempt to build atom group from a leaf atom");
        for(int i=0; i<4; i++) subFactory.makeAtom(group);
        new Configuration4(((LatticeCubicFcc)group).latticeConstant).initializePositions(group);
    }
    
    public double latticeConstant() {return latticeConstant;}
    
    private class Configuration4 extends Configuration {

        public void initializePositions(AtomIterator[] iterators) {
            AtomIterator iterator = iterators[0];
            int i = 0;
            Space.Vector[] dr = unitCell(latticeConstant);
            while(iterator.hasNext()) {
                iterator.next().coord.translateTo(dr[i++]);
            }
        }
    
        private Space.Vector[] unitCell(double latticeConstant){
            Space3D.Vector[] p = new Space3D.Vector[4];
            for(int i=0; i<4; i++) {
                p[i] = new Space3D.Vector();
            }
            p[0].x = 0.0;
            p[0].y = 0.0;
            p[0].z = 0.0;
            p[1].x = 0.0;
            p[1].y = 0.5*latticeConstant;
            p[1].z = 0.5*latticeConstant;
            p[2].x = 0.5*latticeConstant;
            p[2].y = 0.0;
            p[2].z = 0.5*latticeConstant;
            p[3].x = 0.5*latticeConstant;
            p[3].y = 0.5*latticeConstant;
            p[3].z = 0.0;
            return p;
        }//end of unitCell  
    }//end of Configuration4
        
 /*     public static class Criterion1 implements SiteIterator.Neighbor.Criterion{
                Space3D.BoundaryPeriodicSquare periodicBoundary = new Space3D.BoundaryPeriodicSquare(8.0,8.0,8.0);
                public boolean areNeighbors(Site site1,Site site2){
                    double r2 = Space3D.r2(
                        ((BravaisLattice.Coordinate)site1.coordinate()).position(),
                        ((BravaisLattice.Coordinate)site2.coordinate()).position(),
                        periodicBoundary);
                        if(r2 <= 0.5*latticeConstant){
                            return true;
                        }
                        else
                            {return false;}
                            
                }
      }
   */ 
    
}
