package etomica.lattice;
import etomica.*;
import java.util.Observer;
import java.util.Observable;

/**
 * Arbitrary-dimension Bravais Lattice. 
 */
public class BravaisLattice extends Atom implements AbstractLattice {

   Space.Vector[] primitiveVector;
   private double[] primitiveVectorLength;
   private int[] idx;
   private AtomList siteList;
   int D;
   
   private BravaisLattice(Space space, AtomType type) {
        super(space, type);
        D = space.D();
        idx = new int[D];
        primitiveVectorLength = new double[D];
   }
   
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
        for(int i=0; i<D; i++) {
            primitiveVectorLength[i] = Math.sqrt(pVectors[i].squared());
        }
        updateCoordinates();
        notifyObservers();
    }
    public void setPrimitiveVector(int i, Space.Vector pVector) {
        if(pVector == null) return;
        if(i < 0 || i >= D) return;
        primitiveVector[i] = pVector;
        primitiveVectorLength[i] = Math.sqrt(pVector.squared());
        updateCoordinates();
        notifyObservers();
    }
    public void setPrimitiveVector(int vectorIndex, int xyzIndex, double value) {
        primitiveVector[vectorIndex].setComponent(xyzIndex,value);
        primitiveVectorLength[vectorIndex] = Math.sqrt(primitiveVector[vectorIndex].squared());
        updateCoordinates();
        notifyObservers();
    }
    
    //not carefully implemented
    public Site nearestSite(Space.Vector r) {
        for(int i=D-1; i>=0; i--) {
            idx[i] = (int)(r.component(i)/primitiveVectorLength[i] + 0.5);
        }
        return site(idx);
    }
    
    public Site site(int[] idx) {
        Atom site = null;
        int i=0;
        do site = ((AtomTreeNodeGroup)site.node).childList.get(idx[i++]);
        while(!site.node.isLeaf() && i<idx.length);
        return (Site)site;
    }
    
    public AtomList siteList() {return siteList;}
    
    public int D() {return D;}
    
    /**
     * Causes all coordinates to update their position vectors.  
     * Typically invoked when a change is made to the primitive vectors.
     */
    public void updateCoordinates() {
    }

 /**
  * Factory to construct an arbitrary-dimension Bravais lattice.
  * The lattice occupies a D'-dimensional space.  Normally D' = D, but this need not be 
  * the case (e.g., we may have a 2-dimensional (planar) lattice given at some orientation in 3-d).
  * The placement of the sites in the D'-dimensional space is given by the primitiveVector array.
  */
public static class Factory extends AtomFactoryTree {
    
     //dimension of lattice equals the number of primitive vectors
     //dimension of embedding space equals the length of each primitive vector
    public Factory(Space space, int[] dimensions, Space.Vector[] pVectors, AtomFactory siteFactory) {
        super(space, siteFactory, dimensions, configArray(space, pVectors));
        primitiveVectors = pVectors;
    }
    
    public Atom build() {
        BravaisLattice group = new BravaisLattice(space, groupType);
        build(group);
        AtomIteratorTree leafIterator = new AtomIteratorTree(group);
        leafIterator.reset();
        group.siteList = new AtomList(leafIterator);
        return group;
    }

    private static Configuration[] configArray(Space space, Space.Vector[] pVectors) {
        Configuration[] array = new Configuration[pVectors.length];
        for(int i=0; i<array.length; i++) {
            array[i] = new ConfigurationLinear(space);
            ((ConfigurationLinear)array[i]).setOffset(pVectors[i]);
        }
        return array;
    }
    
    private Space.Vector[] primitiveVectors;
}//end of Factory

    /**
     * Main method to demonstrate use of BravaisLattice and to aid debugging
     */
/*    public static void main(String[] args) {
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
    }//end of main
    */
}//end of BravaisLattice