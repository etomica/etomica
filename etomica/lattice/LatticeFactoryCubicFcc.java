package etomica.lattice;
import etomica.*;

public class LatticeFactoryCubicFcc extends AtomFactory {
    
    private double latticeConstant;
    private LatticeFactoryCubic subFactory;
    
    /**
     * Creates an fcc lattice with the given number of sites in a box of the given size,
     * using a simple site factory (Site.Factory).
     */
    public LatticeFactoryCubicFcc(Simulation sim, int nSites, double boxSize) {
        this(sim, getSize(nSites), getLatticeConstant(nSites,boxSize), new Site.Factory()); 
    }
    /**
     * Creates an fcc lattice with the given number of sites in a box of the given size,
     * using the given site factory.
     */
    public LatticeFactoryCubicFcc(Simulation sim, AtomFactory siteFactory, int nSites, double boxSize) {
        this(sim, getSize(nSites), getLatticeConstant(nSites,boxSize), factory); 
    }
    /**
     * Creates an fcc lattice with 4*dimensions[1]*dimensions[2]*dimensions[3] sites, using the given lattice constant
     * and a simple site factory.
     */
    public LatticeFactoryCubicFcc(Simulation sim, int[] dimensions, double latticeConstant) {
        this(sim, new Site.Factory(), dimensions, latticeConstant);
    }
    /**
     * Creates an fcc lattice with 4*dimensions[1]*dimensions[2]*dimensions[3] sites, 
     * using the given lattice constant and site factory.
     */
    public LatticeFactoryCubicFcc(Simulation sim, AtomFactory siteFactory, int[] dimensions, double latticeConstant) {
        super(sim);
        this.latticeConstant = latticeConstant;
        subFactory = new LatticeFactoryCubic(sim, siteFactory, dimensions, latticeConstant);
    }
    
    public Atom build(AtomTreeNodeGroup parent) {
        AtomGroup group = new AtomGroup(parentSimulation.space, groupType, parent);
        return build(group);
    }
    
    public Atom build(Atom group) {
        if(!(group.node instanceof AtomTreeNodeGroup)) throw new IllegalArgumentException("Error in LatticeFactoryCubicFcc.build(Atom):  Attempt to build atom group from a leaf atom");
        for(int i=0; i<4; i++) subFactory.makeAtom(group);
        configuration.initializePositions(group);
    }
    
    public double latticeConstant() {return latticeConstant;}
    
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


    
    private static Space.Vector[] unitCell(double latticeConstant){
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
