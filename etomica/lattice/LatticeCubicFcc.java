package etomica.lattice;
import etomica.Atom;
import etomica.AtomFactory;
import etomica.AtomIterator;
import etomica.AtomIteratorList;
import etomica.AtomIteratorTree;
import etomica.AtomList;
import etomica.AtomSequencer;
import etomica.AtomSequencerFactory;
import etomica.AtomTreeNodeGroup;
import etomica.AtomType;
import etomica.Configuration;
import etomica.Simulation;
import etomica.SimulationEventManager;
import etomica.Space;
import etomica.Space3D;
import etomica.action.AtomAction;
import etomica.action.AtomActionAdapter;

/**
 * fcc lattice formed to a cubic shape.  Note: does not extend BravaisLattice.
 */

/* History of changes
 * 09/18/02 (DAK) modified site method to return Atom instead of Site.
 */
public class LatticeCubicFcc extends Atom implements AbstractLattice {
    
   private double latticeConstant;
    
   private final AtomList siteList = new AtomList();
   private int[] dimensions;
   private final int D = 3;
   private int[] idx = new int[D];
   private NeighborManager.Criterion neighborCriterion;
   private final SimulationEventManager eventManager = new SimulationEventManager();
   private final LatticeEvent rebuildEvent = new LatticeEvent(this, LatticeEvent.REBUILD);
   private final LatticeEvent allSiteEvent = new LatticeEvent(this, LatticeEvent.ALL_SITE);
   private final LatticeEvent resetNbrEvent = new LatticeEvent(this, LatticeEvent.RESET_NBRS);
   
   /**
    * Constructor should be invoked only within a build() method of the
    * Factory class.  The build method handles the construction of the
    * substructure under this instance, which forms the lattice.
    */
   private LatticeCubicFcc(Space space, AtomType type, int[] dimensions, AtomTreeNodeGroup parent,
                            double a) {
        super(space, type, AtomTreeNodeGroup.FACTORY, 
                AtomSequencerFactory.SIMPLE, parent);
        if(space.D() != 3) throw new IllegalArgumentException("Error in LatticeCubicFcc constructor:  Given space is not 3D");
        this.dimensions = new int[dimensions.length];
        latticeConstant = a;
        System.arraycopy(dimensions, 0, this.dimensions, 0, dimensions.length);
   }
   
    /**
     * Creates an fcc lattice with the given number of sites in a box of the given size,
     * using a simple site factory (Site.Factory).
     */
    public static LatticeCubicFcc makeLattice(Space space, int nSites, double boxSize) {
        return makeLattice(space, new Site.Factory(space), calculateSize(nSites), 
                            calculateLatticeConstant(nSites,boxSize)); 
    }
    /**
     * Creates an fcc lattice with the given number of sites in a box of the given size,
     * using the given site factory.
     */
    public static LatticeCubicFcc makeLattice(Space space, AtomFactory siteFactory, 
                                                int nSites, double boxSize) {
        return makeLattice(space, siteFactory, calculateSize(nSites), 
                            calculateLatticeConstant(nSites,boxSize)); 
    }
    /**
     * Creates an fcc lattice with 4*dimensions[1]*dimensions[2]*dimensions[3] sites, using the given lattice constant
     * and a simple site factory.
     */
    public static LatticeCubicFcc makeLattice(Space space, int[] dimensions, double latticeConstant) {
        return makeLattice(space, new Site.Factory(space), dimensions, latticeConstant);
    }
    
    public static LatticeCubicFcc makeLattice(Space space, AtomFactory siteFactory, 
                                                int[] dimensions, double latticeConstant) {
        return (LatticeCubicFcc)new Factory(space, AtomSequencerFactory.SIMPLE, siteFactory, dimensions, latticeConstant).makeAtom();
    }
    
    /**
     * Determines the size of the cubic lattice consistent with an fcc lattice having
     * the given number of sites.  
     * Returns the smallest k such that 4*k^3 is greater than or equal to n.
     */
    private static int[] calculateSize(int n) {
        if(n <= 0) return new int[] {0, 0, 0};
        int nLat = 0;
        int k = 0;
        while(n > nLat) {
            k++;
            nLat = 4*k*k*k;
        }
        return new int[] {k, k, k};
    }
    
    private static double calculateLatticeConstant(int nSites, double boxSize) {
        int n = calculateSize(nSites)[0];
        return boxSize/(double)n;
    }

    /**
     * Returns the event manager the registers listeners and notifies them
     * of events indicating changes in this lattice.
     */
    public SimulationEventManager eventManager() {return eventManager;}

    public int[] getDimensions() {return dimensions;}
    public void setDimensions(int[] dim) {
        if(dim.length != 3) throw new IllegalArgumentException("Incorrect length of dimensions array given to LatticeCubicFcc.setDimensions");
        for(int i=0; i<dim.length; i++) {
            dimensions[i] = dim[i];
        }
        rebuild();
    }
    
    public void setLatticeConstant(double a) {
        if(a < 0) throw new IllegalArgumentException("Invalid (negative) value given to LatticeCubicFcc.setLatticeConstant");
        latticeConstant = a;
        update();
    }
    public double getLatticeConstant() {return latticeConstant;}
    
    
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
     * Returns the site corresponding to the given index.  First index entry
     * indicates one of the four simple-cubic sublattices (acceptable values
     * for this entry are 0, 1, 2, or 3); the remaining three entries indicate
     * the element of selected simple-cubic sublattice. 
     */
     //needs work to improve probable inefficiency with list.get
     //and to handle index values out of size of lattice
    public Atom site(int[] idx) {
        if(idx.length != 4) throw new IllegalArgumentException("Length of array given to LatticeCubicFcc.site(int[]) is incorrect; must be of length 4");
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
    
    public int siteCount() {return siteList.size();}
    
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
     * Reconstructs all sites of the lattice.  Invoked when dimensions of
     * lattice are changed.
     */
    public void rebuild() {
        creator().build(this);
        eventManager.fireEvent(rebuildEvent);
    }
     
    
    /**
     * Causes all coordinates to update their position vectors, and
     * notifies any observers that a change has occurred.
     */
    public void update() {
        creator().getConfiguration().initializePositions(this);
        eventManager.fireEvent(allSiteEvent);
//        notifyObservers();
    }//end of update
    
    /**
     * Sets up neighbors using the default criterion, which is the criterion given
     * in the most recent call to setupNeighbors(Criterion).  If this was not previously
     * invoked, the adjacency criterion is used.
     */
    public void setupNeighbors() {
        if(neighborCriterion == null) neighborCriterion = new AdjacencyCriterion(this);
        setupNeighbors(neighborCriterion);
    }
    public void setupNeighbors(NeighborManager.Criterion criterion) {
        AtomIteratorList iterator = new AtomIteratorList(siteList);
        iterator.reset();
        while(iterator.hasNext()) {
            Site site = (Site)iterator.nextAtom();
            site.neighborManager().setupNeighbors(siteList, criterion);
        }
        eventManager.fireEvent(resetNbrEvent);
    }//end of setupNeighbors
    
    public Space.Vector[] positions() {
        int n = siteCount();
        Space.Vector[] r = new Space.Vector[n];
        for(int i=0; i<n; i++) {r[i] = new Space3D.Vector();}
        AtomIteratorList iteratorsites = new AtomIteratorList(siteList);
        iteratorsites.reset();
        int i = 0;
        while (iteratorsites.hasNext()){
            Atom site = iteratorsites.nextAtom();
            r[i++].E(site.coord.position());
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
    private LatticeCubicFcc lattice;
    
    public AdjacencyCriterion(LatticeCubicFcc lattice) {
        this.lattice = lattice;
        
        //need to rewrite areNeighbors method
        throw new RuntimeException("LatticeCubicFcc.AdjacencyCriterion not yet implemented");
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
        int[] dim = lattice.getDimensions();
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
    private double latticeConstant;
    private int[] dimensions;
    private AtomFactory siteFactory;
    
    /**
     * Creates an fcc lattice with 4*dimensions[1]*dimensions[2]*dimensions[3] sites, 
     * using the given lattice constant and site factory.
     */
    public Factory(Space space, AtomSequencerFactory seqFactory, AtomFactory siteFactory, int[] dimensions, double latticeConstant) {
        super(space, seqFactory);
        this.siteFactory = siteFactory;
        this.dimensions = new int[dimensions.length];
        System.arraycopy(dimensions, 0, this.dimensions, 0, dimensions.length);
        this.latticeConstant = latticeConstant;
        configuration = new Configuration4(space,latticeConstant);
    }
    
    public boolean isGroupFactory() {return false;}
    
    /**
     * Instantiates and builds a new fcc lattice using the latticeConstant, dimensions,
     * and atomFactory as currently set in this factory.
     */
    public Atom build(AtomTreeNodeGroup parent) {
        groupType.childrenAreGroups = true;
        Atom group = new LatticeCubicFcc(space, groupType, dimensions, parent, latticeConstant);
        return build(group);
    }
    
    /**
     * Builds the given atom into an fcc lattice.  Uses the current values of latticeConstant
     * and dimensions of the atom, and constructs sites using this factory's atomFactory.
     */
    public Atom build(Atom atom) {
        if(!(atom instanceof LatticeCubicFcc)) throw new IllegalArgumentException("Error in LatticeFactoryCubicFcc.build(Atom):  Attempt to build atom group from a leaf atom");
        LatticeCubicFcc group = (LatticeCubicFcc)atom;
        ((AtomTreeNodeGroup)group.node).removeAllChildren();
        LatticeFactoryCubic subFactory = 
            new LatticeFactoryCubic(space, siteFactory, group.getDimensions(), group.getLatticeConstant());
        for(int i=0; i<4; i++) subFactory.makeAtom((AtomTreeNodeGroup)group.node);
        //set up siteList
        AtomIteratorTree leafIterator = new AtomIteratorTree();
        leafIterator.setRoot(group);
        leafIterator.reset();
        group.siteList.clear();
        group.siteList.addAll(leafIterator);
        //position sites
        configuration.initializePositions(group);
        return atom;
    }
    
    public double getLatticeConstant() {return latticeConstant;}
    public void setLatticeConstant(double a) {latticeConstant = a;}
    
    private static class Configuration4 extends Configuration {

        private double latticeConstant;
        Configuration4(Space space, double a) {
            super(space);
            latticeConstant = a;
        }
        
        public void initializePositions(Atom atom) {
            LatticeCubicFcc lattice = (LatticeCubicFcc)atom;
            latticeConstant = lattice.getLatticeConstant();
            super.initializePositions(atom);
        }
        
        public void initializePositions(AtomIterator[] iterators) {
            AtomIterator iterator = iterators[0];
            int i = 0;
            Space.Vector[] dr = unitCell(latticeConstant);
            while(iterator.hasNext()) {
                iterator.nextAtom().coord.translateTo(dr[i++]);
            }
        }
    
        private Space.Vector[] unitCell(double latticeConstant){
            Space3D.Vector[] p = new Space3D.Vector[4];
            for(int i=0; i<4; i++) {
                p[i] = new Space3D.Vector();
            }
            p[0].setX(0,0.0);
            p[0].setX(1,0.0);
            p[0].setX(2,0.0);
            
            p[1].setX(0,0.0);
            p[1].setX(1,0.5*latticeConstant);
            p[1].setX(2,0.5*latticeConstant);
            
            p[2].setX(0,0.5*latticeConstant);
            p[2].setX(1,0.0);
            p[2].setX(2,0.5*latticeConstant);
            
            p[3].setX(0,0.5*latticeConstant);
            p[3].setX(1,0.5*latticeConstant);
            p[3].setX(2,0.0);
            return p;
        }//end of unitCell  
    }//end of Configuration4
}//end of Factory   
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
   
    /**
     * Main method to demonstrate use of BravaisLattice and to aid debugging
     */
    public static void main(String[] args) {
        System.out.println("main method for LatticeCubicFcc");
        Space space = new Space3D();
        Simulation sim = new Simulation(space);
        Simulation.instance = sim;
        int D = space.D();
        final int nx = 2, ny = 1, nz = 2;
        LatticeCubicFcc lattice = LatticeCubicFcc.makeLattice(space, 
                                new Site.Factory(space),
                                new int[] {nx,ny,nz},
                                1.0);        
        System.out.println("Total number of sites: "+lattice.siteList().size());
        System.out.println();
        System.out.println("Coordinate printout");
        AtomIteratorList iterator = new AtomIteratorList(lattice.siteList());
        iterator.reset();
        while(iterator.hasNext()) {  //print out coordinates of each site
            System.out.print(iterator.nextAtom().coord.position().toString()+" ");
        }
        System.out.println();
        
        System.out.println("Same, using allAtoms method");
        AtomAction printSites = new AtomActionAdapter() {
            public void actionPerformed(Atom s) {
                System.out.print(s.coord.position().toString()+" ");
       //         System.out.println(((Site)s).latticeCoordinate()[1]);
                if(((Site)s).latticeCoordinate()[1]==ny-1) System.out.println();
            }
        };
 //       iterator.allAtoms(printSites);
        System.out.println();
        
        System.out.println();
        System.out.println("Changing dimensions");
        lattice.setDimensions(new int[] {nx, ny+1, nz});
//        iterator.allAtoms(printSites);
        System.out.println();
        
        Atom testSite = lattice.site(new int[] {0,1,0,1});
        int[] idx = null;
        
        System.out.println();
        System.out.println("Translating lattice by (-1.0, 2.0, 0.0)");
        lattice.coord.translateBy(Space.makeVector(new double[] {-1.0, 2.0, 0.0}));
//        iterator.allAtoms(printSites);
        System.out.println();
 
        System.out.println();
        System.out.println("Translating to origin");
        lattice.shiftFirstToOrigin();
//        iterator.allAtoms(printSites);
        System.out.println();
   /*     
        System.out.print("Accessing site (1,1): ");
        testSite = lattice.site(new int[] {1,1});
        idx = ((Site)testSite).latticeCoordinate();
        System.out.println(testSite.toString()+testSite.coord.position().toString());
        System.out.print("latticeCoordinate: ");
        for(int i=0; i<idx.length; i++) System.out.print(idx[i]);
        System.out.println();
        System.out.println();
        
        System.out.print("Accessing site (2,0): ");
        testSite = lattice.site(new int[] {2,0});
        idx = ((Site)testSite).latticeCoordinate();
        System.out.println(testSite.toString()+testSite.coord.position().toString());
        System.out.print("latticeCoordinate: ");
        for(int i=0; i<idx.length; i++) System.out.print(idx[i]);
        System.out.println();
        System.out.println();
 /*       
        lattice.setupNeighbors(new AdjacencyCriterion(lattice));
        SiteIteratorNeighbor nbrIterator = new SiteIteratorNeighbor();
        nbrIterator.setBasis(testSite);

        System.out.println("Sites up-neighbor to this site:");
        nbrIterator.reset(IteratorDirective.UP);
        while(nbrIterator.hasNext()) {  //print out coordinates of each site
            System.out.print(nbrIterator.next().toString()+" ");
        }
        System.out.println();
        System.out.println();
        
        System.out.println("Sites down-neighbor to this site:");
        nbrIterator.reset(IteratorDirective.DOWN);
        while(nbrIterator.hasNext()) {  //print out coordinates of each site
            System.out.print(nbrIterator.next().toString()+" ");
        }
        System.out.println();
        System.out.println();
        
        System.out.println("All neighbors of this site:");
        nbrIterator.reset(IteratorDirective.BOTH);
        while(nbrIterator.hasNext()) {  //print out coordinates of each site
            System.out.print(nbrIterator.next().toString()+" ");
        }
        System.out.println();
        System.out.println();

        System.out.print("A randomly selected site: ");
        testSite = lattice.siteList().getRandom();
        idx = ((Site)testSite).latticeCoordinate();       
        System.out.println(testSite.toString());
        System.out.print("latticeCoordinate: ");
        for(int i=0; i<idx.length; i++) System.out.print(idx[i]);
        System.out.println();
        
        nbrIterator.setBasis(testSite);

        System.out.println("Sites up-neighbor to this site:");
        nbrIterator.reset(IteratorDirective.UP);
        while(nbrIterator.hasNext()) {  //print out coordinates of each site
            System.out.print(nbrIterator.next().toString()+" ");
        }
        System.out.println();
        System.out.println();
        
        System.out.println("Sites down-neighbor to this site:");
        nbrIterator.reset(IteratorDirective.DOWN);
        while(nbrIterator.hasNext()) {  //print out coordinates of each site
            System.out.print(nbrIterator.next().toString()+" ");
        }
        System.out.println();
        System.out.println();
        
        System.out.println("All neighbors of this site:");
        nbrIterator.reset(IteratorDirective.BOTH);
        while(nbrIterator.hasNext()) {  //print out coordinates of each site
            System.out.print(nbrIterator.next().toString()+" ");
        }
        System.out.println();
        System.out.println();
 */      
    }//end of main
   
    
}//end of LatticeCubicFcc
