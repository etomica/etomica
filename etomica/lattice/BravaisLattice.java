package etomica.lattice;
import etomica.*;
import java.util.Observer;
import java.util.Observable;

/**
 * Arbitrary-dimension Bravais Lattice. 
 */
public class BravaisLattice extends Atom implements AbstractLattice {

   private Space.Vector[] primitiveVector;
   private double[] primitiveVectorLength;
   private int[] idx;
   private AtomList siteList;
   int D;
   
   protected BravaisLattice(Space space, AtomType type) {
        super(space, type);
        D = space.D();
        idx = new int[D];
        primitiveVectorLength = new double[D];
        primitiveVector = new Space.Vector[D];
        for(int i=0; i<D; i++) {primitiveVector[i] = space.makeVector();}
   }
   
   /**
    * Constructs a unique BravaisLattice factory and returns a new lattice from it.
    */
   public static BravaisLattice makeLattice(Space space, int[] dimensions, Space.Vector[] pVectors, AtomFactory siteFactory) {
        return (BravaisLattice)new Factory(space, dimensions, pVectors, siteFactory).build();
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
        if(pVectors.length != D) throw new IllegalArgumentException("Error in BravaisLattice.setPrimitiveVector:  incorrect length for primitive-vector array");
        Atom a = this;
        for(int i=0; i<D; i++) {
           primitiveVector[i].E(pVectors[i]);
           ((ConfigurationLinear)a.creator().getConfiguration()).setOffset(pVectors[i]);
           a = a.node.firstChildAtom();
        }
        type.creator().getConfiguration().initializePositions(this);
        notifyObservers();
    }
    public void setPrimitiveVector(int i, Space.Vector pVector) {
        if(pVector == null) return;
        ((ConfigurationLinear)getAtomType(i).creator().getConfiguration()).setOffset(pVector);
        primitiveVector[i].E(pVector);
        type.creator().getConfiguration().initializePositions(this);
        notifyObservers();
    }
    
    /**
     * Returns the atom type instance held by all atoms at the depth d.
     */
    protected AtomType getAtomType(int d) {
        if(d < 0 || d >= D) throw new IllegalArgumentException("Error in BravaisLattice.getFactory: inappropriate dimension specified");
        Atom a = this;
        int i = d;
        while(i < d) {a = a.node.firstChildAtom(); i++;}
        return a.type;
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
     //copies of primitive vectors are made during construction of lattice, so subsequent alteration
     //of them by the calling program has no effect on lattice vectors
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
        group.setPrimitiveVector(primitiveVectors);
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
    public static void main(String[] args) {
        System.out.println("main method for BravaisLattice");
        Space space = new Space2D();
        int D = space.D();
        BravaisLattice lattice = BravaisLattice.makeLattice(space, 
                                new int[] {3,2},
                                new Space.Vector[] {Space.makeVector(new double[] {1.,0.}),
                                                   Space.makeVector(new double[] {0.,1.}),},
                                new AtomFactoryMono(space));        
        System.out.println("Total number of sites: "+lattice.siteList().size());
        System.out.println();
        System.out.println("Coordinate printout");
        AtomIteratorList iterator = new AtomIteratorList(lattice.siteList());
        iterator.reset();
        while(iterator.hasNext()) {  //print out coordinates of each site
            System.out.print(iterator.next().toString()+" ");
        }
        System.out.println();
        
        System.out.println("Same, using allAtoms method");
        AtomAction printSites = new AtomAction() {public void actionPerformed(Atom s) {System.out.print(s.toString()+" ");}};
        iterator.allAtoms(printSites);
        System.out.println();
        
        System.out.println();
        System.out.println("Changing primitive vector");
        lattice.setPrimitiveVector(new Space.Vector[] {Space.makeVector(new double[] {0.,1.}),
                                                       Space.makeVector(new double[] {0.5,0.})});
        iterator.allAtoms(printSites);
        System.out.println();
        
        System.out.println();
        System.out.println("Translating lattice by (-1.0, 2.0)");
        lattice.coord.translateBy(Space.makeVector(new double[] {-1.0, 2.0}));
        iterator.allAtoms(printSites);
        System.out.println();
        
        System.out.print("Accessing site (1,1): ");
        Site testSite = lattice.site(new int[] {1,1});
        System.out.println(testSite.toString());
        System.out.println();
    /*            
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
*/
        System.out.print("A randomly selected site: ");
        System.out.println(lattice.siteList().getRandom().toString());
    }//end of main
}//end of BravaisLattice