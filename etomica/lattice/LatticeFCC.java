package etomica.lattice;
import etomica.Space3D;

public class LatticeFCC extends LatticeCubic {
    
    /**
     * Creates an fcc lattice with the given number of sites in a box of the given size,
     * using a simple site factory (Site.Factory).
     */
    public LatticeFCC(int nSites, double boxSize) {
        this(getSize(nSites), getLatticeConstant(nSites,boxSize), new Site.Factory()); 
    }
    /**
     * Creates an fcc lattice with the given number of sites in a box of the given size,
     * using the given site factory.
     */
    public LatticeFCC(int nSites, double boxSize, SiteFactory factory) {
        this(getSize(nSites), getLatticeConstant(nSites,boxSize), factory); 
    }
    /**
     * Creates an fcc lattice with 4*size[1]*size[2]*size[3] sites, using the given lattice constant
     * and a simple site factory.
     */
    public LatticeFCC(int[] size, double latticeConstant) {
        this(size, latticeConstant, new Site.Factory());
    }
    /**
     * Creates an fcc lattice with 4*size^3 sites, using the given lattice constant
     * and site factory.
     */
    public LatticeFCC(int[] size, double latticeConstant, SiteFactory factory) {
        super(size, latticeConstant, new Basis(positions(latticeConstant),factory));
    }
    
/*    public LatticeFCC(){
        super(3, 4, primaryVectors(), basis1());
        //set up neighbors
        SiteIterator.Neighbor iterator = new SiteIterator.Neighbor();
        iterator.reset();
        SiteIterator.Neighbor.Cursor cursor =  iterator.makeCursor();
        Criterion1 criterion1 = new Criterion1() ;
        while(iterator.hasNext()) {
            cursor.reset();
            Site site = iterator.next();
            site.adjacentIterator().setNeighbors(cursor, criterion1);
        }
    }
 */   
    private static Space3D.Vector[] positions(double latticeConstant){
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