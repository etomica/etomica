package etomica.lattice;
import etomica.Space;
import java.util.Observer;
import java.util.Observable;
import java.util.Arrays;

/**
 * Lattice with a basis
 */
public class Lattice implements AbstractLattice, Observer {
    
    private BravaisLattice lattice;
    private Basis basis;
    private SiteIterator.List iterator = new SiteIterator.List();
    
    /**
     * Create a square Bravais lattice (same number of sites in each dimension) of dimension "d" with "size" elements in each dimension
     * Need java 1.2 for Array.fill
     */
    public Lattice(int d, int size, Space.Vector[] pVectors, Basis basis) {
        this(filledArray(d,size), pVectors, basis);
    }
    /**
     * Creates a Bravais lattice with the given basis.  Number of lattice sites in each
     * dimension is determined by the given integer array.  Size of this array also
     * determines the dimension of the lattice.
     */
    public Lattice(int[] dimensions, Space.Vector[] pVectors, Basis basis) {
        lattice = new BravaisLattice(dimensions, basis, pVectors);
        this.basis = basis;
        lattice.addObserver(this); //to update site positions if primitive vectors change
        basis.addObserver(this);   //to update site positions if basis vectors change
        //Construct iterator of all sites
        SiteIterator lIterator = lattice.iterator();
        lIterator.reset();
        while(lIterator.hasNext()) { //loop over all basis points
            Basis.SiteLattice b = (Basis.SiteLattice)lIterator.next();
            SiteIterator bIterator = b.iterator();
            bIterator.reset();
            while(bIterator.hasNext()) { //collect sites from each basis
                iterator.addSite(bIterator.next());
            }
        }
        updateCoordinates();
    }
    
    /**
     * Updates position coordinates of all sites on the lattice.
     * Positions are computed using current values of primitive and basis vectors.
     * The positions must be updated if these vectors are changed.  Updating is performed
     * automatically by the accessor methods of the vectors in Basis and BravaisLattice.
     */
    public void updateCoordinates() {
        iterator.reset();
        while(iterator.hasNext()) {
            ((Basis.SiteLattice.Coordinate)iterator.next().coordinate()).update();
        }
    }
    /**
     * Implementation of Observer interface; updates position coordinates of all sites in lattice.
     * Notification of need for update comes from lattice or basis, if their defining vectors are changed.
     */
    public void update(Observable o, Object arg) {updateCoordinates();}
 
 //  construct array of int with d elements each having value size;
   private static int[] filledArray(int d, int size) {
     int[] array = new int[d];
     Arrays.fill(array, size);
     return array;
   }
   
  //Lattice interface methods
    public SiteIterator iterator() {iterator.reset(); return iterator;}
    public int D() {return lattice.D();}           
    public int siteCount() {return lattice.siteCount() * basis.siteCount();} 
//    public int coordinationNumber() {return lattice.coordinationNumber();}   

//not completed
    public Site site(AbstractLattice.Coordinate coord) {
//        Coordinate coordinate = (Coordinate)coord;
//        Basis.SiteLattice basisSite = (Basis.SiteLattice)lattice.site(coordinate.latticeCoordinate);
//        return basisSite.site(coordinate.basisCoordinate);
            return null;
    }
    public Site randomSite() {return ((Basis.SiteLattice)lattice.randomSite()).randomSite();}            
//    public Site.Iterator iterator() {return lattice.iterator();} 


    public BravaisLattice lattice() {return lattice;}
    public Basis basis() {return basis;}
    
 /*   public class Coordinate implements AbstractLattice.Coordinate {
        private AbstractLattice.Coordinate latticeCoordinate;
        private AbstractLattice.Coordinate basisCoordinate;
        public Coordinate(AbstractLattice.Coordinate lCoord, AbstractLattice.Coordinate bCoord) {
            latticeCoordinate = lCoord;
            basisCoordinate = bCoord;
        }
    }*/
}