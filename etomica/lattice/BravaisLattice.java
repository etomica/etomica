package etomica.lattice;
import etomica.Space;
import java.util.Random;
import java.util.Observer;
import java.util.Observable;

/**
 * Arbitrary-dimension Bravais Lattice.  
 * The lattice occupies a D'-dimensional space.  Normally D' = D at the top level, but this need not be 
 * the case (e.g., we may have a 2-dimensional (planar) lattice given at some orientation in 3-d).
 * The placement of the sites in the D'-dimensional space is given by the primitiveVector array.  
 */
public class BravaisLattice extends IntegerBravaisLattice { ///implements AbstractLattice {

   Space.Vector[] primitiveVector;
   Space.Vector translation;  //amount that lattice has been displaced via call to translateBy
   //anonymous inner class passed to superclass IntegerBravaisLattice constructor to
   //have it use the BravaisLattice coordinates for placement in the sites.
   private static final IntegerBravaisLattice.IndexedCoordinate.Factory coordFactory = 
        new IntegerBravaisLattice.IndexedCoordinate.Factory() {
                public IndexedCoordinate makeCoordinate(int[] i) {return new Coordinate(i);}
             };
   
     //dimension of lattice equals the number of primitive vectors
     //dimension of embedding space equals the length of each primitive vector
   public BravaisLattice(int[] dimensions, SiteFactory factory, Space.Vector[] pVectors) {
        super(dimensions, factory, coordFactory);
   ///     super(dimensions, factory, coordFactory(dimensions.length));
        primitiveVector = pVectors;
        translation = Space.makeVector(D());
        //Set lattice for each site
        SiteIterator iter = iterator();
        iter.reset();
        while(iter.hasNext()) {
            Site site = iter.next();
            ((Coordinate)site.coordinate()).setLattice(this);
        }
   }
   
///   private static IntegerBravaisLattice.IndexedCoordinate.Factory coordFactory(final int D) {
///        return new IntegerBravaisLattice.IndexedCoordinate.Factory() {
///                public IndexedCoordinate makeCoordinate(int[] i) {return new Coordinate(i,D);}
///             };
///   }
   
    /**
     * Observable anonymous inner class that notifies interested objects if the primitive vectors change.
     * The observer would then update the sites created with this basis.  Used by the Lattice class.
     */
    private Observable monitor = new Observable() {
        public void notifyObservers() {this.notifyObservers(null);}
        public void notifyObservers(Object obj) {
            setChanged();
            super.notifyObservers(obj);
        }
        public void addObserver(Observer o) {if(o != null) super.addObserver(o);}
    };
    /**
     * Registers the given Observer so that it is notified any time the primitive vectors are modified.
     */
    public void addObserver(Observer o) {monitor.addObserver(o);}
    /**
     * Notifies observers that the primitive vectors are modified.
     */
    protected void notifyObservers() {monitor.notifyObservers();}

    public void setPrimitiveVector(Space.Vector[] pVectors) {
        if(pVectors == null) return;
        primitiveVector = pVectors;
        updateCoordinates();
        notifyObservers();
    }
    public void setPrimitiveVector(int i, Space.Vector pVector) {
        if(pVector == null) return;
        if(i < 0 || i >= D) return;
        primitiveVector[i] = pVector;
        updateCoordinates();
        notifyObservers();
    }
    public void setPrimitiveVector(int vectorIndex, int xyzIndex, double value) {
        primitiveVector[vectorIndex].setComponent(xyzIndex,value);
        updateCoordinates();
        notifyObservers();
    }
    /**
     * Translates all the sites in the lattice by the given vector
     */
    public void translateBy(Space.Vector d) {
        SiteIterator iterator = iterator();
        iterator.reset();
        while(iterator.hasNext()) {
            ((Coordinate)iterator.next().coordinate()).position().PE(d);
        }
        translation.PE(d);
    }
    
    /**
     * Causes all coordinates to update their position vectors.  
     * Typically invoked when a change is made to the primitive vectors.
     */
    public void updateCoordinates() {
        SiteIterator iterator = iterator();
        iterator.reset();
        while(iterator.hasNext()) {
            ((Coordinate)iterator.next().coordinate()).update();
        }
    }
    
    public static class Coordinate extends IntegerBravaisLattice.Coordinate implements AbstractLattice.PositionCoordinate {
        Space.Vector r;
        BravaisLattice bravaisLattice;
        public Coordinate(int[] index) {
            super(index);
        }
 ///       public Coordinate(int[] index, int D) {
 ///           super(index);
 ///           r = Space.makeVector(D);
 ///       }
        public void setLattice(BravaisLattice lattice) {
            if(bravaisLattice != null) return; //can set lattice only once
            bravaisLattice = lattice;
            r = Space.makeVector(bravaisLattice.D());
            update();
        }
        public final void update() {
            r.E(0.0);
            for(int j=0; j<bravaisLattice.D(); j++) { //loop over primitive vectors
                r.PEa1Tv1((double)index()[j],bravaisLattice.primitiveVector[j]);
            }
        }
        public Space.Vector position() {return r;}
        public String toString() {
            String value = "( ";
            for(int i=0; i<r.length(); i++) value += r.component(i) + " ";
            value += ") ";
            return value;
        }
    
    }//end of BravaisLattice.Coordinate
   
    /**
     * Main method to demonstrate use of BravaisLattice and to aid debugging
     */
    public static void main(String[] args) {
        System.out.println("main method for BravaisLattice");
        int D = 2;
        BravaisLattice lattice = 
            new BravaisLattice(new int[] {3,2}, new Site.Factory(),
                               new Space.Vector[] {Space.makeVector(new double[] {1.,0.}),
                                                   Space.makeVector(new double[] {0.,1.})});        
        System.out.println("Total number of sites: "+lattice.siteCount());
        System.out.println();
        System.out.println("Coordinate printout");
        SiteIterator iterator = lattice.iterator();
        iterator.reset();
        while(iterator.hasNext()) {  //print out coordinates of each site
            System.out.print(iterator.next().toString()+" ");
        }
        System.out.println();
        
        System.out.println("Same, using allSites method");
        SiteAction printSites = new SiteAction() {public void actionPerformed(Site s) {System.out.print(s.toString()+" ");}};
        iterator.allSites(printSites);
        System.out.println();
        
        System.out.println();
        System.out.println("Changing primitive vector");
        lattice.setPrimitiveVector(new Space.Vector[] {Space.makeVector(new double[] {0.,1.}),
                                                       Space.makeVector(new double[] {0.5,0.})});
        iterator.allSites(printSites);
        System.out.println();
        
        System.out.println();
        System.out.println("Translating lattice by (-1.0, 2.0)");
        lattice.translateBy(Space.makeVector(new double[] {-1.0, 2.0}));
        iterator.allSites(printSites);
        System.out.println();
        
        System.out.print("Accessing site (1,1): ");
        Site testSite = lattice.site(new IntegerBravaisLattice.Coordinate(new int[] {1,1}));
        System.out.println(testSite.toString());
        System.out.println();
                
        System.out.println("Sites up-neighbor to this site:");
        iterator = testSite.neighborIterator();
        ((SiteIterator.Neighbor)iterator).resetUp();
        while(iterator.hasNext()) {  //print out coordinates of each site
            System.out.print(iterator.next().toString()+" ");
        }
        System.out.println();
        System.out.println();
        
        System.out.println("Sites down-neighbor to this site:");
        iterator = testSite.neighborIterator();
        ((SiteIterator.Neighbor)iterator).resetDown();
        while(iterator.hasNext()) {  //print out coordinates of each site
            System.out.print(iterator.next().toString()+" ");
        }
        System.out.println();
        System.out.println();

        System.out.print("A randomly selected site: ");
        System.out.println(lattice.randomSite().toString());
    }
}