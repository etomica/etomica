package etomica.lattice;
import etomica.*;
import etomica.action.AtomActionAdapter;

/**
 * Arbitrary-dimension Bravais Lattice. 
 */
 
 /* History of changes
  * 09/16/02 (DAK) got rebuild() method working by invoking node.removeAllChildren() 
  *           before sending lattice to factory for rebuilding.
  *           added accessor methods for Factory dimensions field (changes made to
  *           AtomFactoryTree permit this field to be changed.
  * 09/18/02 (DAK) modified site method to return Atom instead of Site.
  * 01/12/03 (DAK) commented/uncommented parts of setupNeighbors(Criterion), as
  * indicated there
  * 01/19/04 (DAK) added leafList method; modified in accord with change
  * involving Basis class
  * 01/20/04 (DAK) revised makeLattice methods and Factory class constructors to
  * take a Crystal as an argument.
  */
public class BravaisLattice extends Atom implements AbstractLattice {

   private Factory factory;
   private double[] primitiveVectorLength;
   private int[] idx;
   private final AtomList siteList = new AtomList();
   private boolean delayUpdate = false;//used to signal update and rebuild methods to not perform update, because another call will be forthcoming
   int D;
   private NeighborManager.Criterion neighborCriterion;
   public final SimulationEventManager eventManager = new SimulationEventManager();
   private final LatticeEvent rebuildEvent = new LatticeEvent(this, LatticeEvent.REBUILD);
   private final LatticeEvent allSiteEvent = new LatticeEvent(this, LatticeEvent.ALL_SITE);
   private final LatticeEvent resetNbrEvent = new LatticeEvent(this, LatticeEvent.RESET_NBRS);
   
   /**
    * Constructor should be invoked only within a build() method of a
    * Factory class.  The build method handles the construction of the
    * tree structure under this instance, which forms the lattice.
    */
   private BravaisLattice(Space space, AtomType type, AtomTreeNodeGroup parent) {
        super(space, type, AtomTreeNodeGroup.FACTORY, 
                AtomSequencerSimple.FACTORY, parent);
        D = space.D();
        idx = new int[D];
        primitiveVectorLength = new double[D];
   }
   
   /**
    * Constructs and returns a BravaisLattice.
    */
    public static BravaisLattice makeLattice(Space space, Crystal crystal, int[] dimensions) {
		return (BravaisLattice)new Factory(space, crystal, dimensions).makeAtom();
    }
    /**
     * Constructs and returns BravaisLattice instance of default size as square
     * of 3 units in each direction.
     * @param space
     * @param dimensions
     * @param primitive
     * @return BravaisLattice
     */
    public static BravaisLattice makeLattice(Space space, Crystal crystal) {
    	return makeLattice(space, crystal, space.makeArrayD(3));
    }
	/**
	 * Constructs and returns a BravaisLattice with a single-atom basis.
	 */
	 public static BravaisLattice makeLattice(Space space, int[] dimensions, Primitive primitive) {
		 return makeLattice(space, new Crystal(primitive, new Basis(space)), dimensions);
	 }
   /**
    * Constructs and returns a BravaisLattice with each site made by the given
    * AtomFactory.
    */
    public static BravaisLattice makeLattice(
                Space space, AtomFactory siteFactory, int[] dimensions, Primitive primitive) {
        return makeLattice(space, new Crystal(primitive, new Basis(space, 1, siteFactory)), dimensions);
    }
               
    public int[] getDimensions() {return factory.dimensions;}
    public void setDimensions(int[] dim) {
        factory.setDimensions(dim);
        rebuild();
    }
    
    public Atom centerSite() {
    	int[] halfDim = new int[D];
    	for(int i=0; i<D; i++) halfDim[i] = (int)factory.dimensions[i]/2;
		return site(halfDim);
    }
    
    /**
     * Sets the size of the lattice (number of atoms in each direction) so that
     * it has a number of sites equal or greater than the given value n.  Sets
     * lattice to be square  (same number in all directions), and finds smallest
     * size that gives number of sites equal or exceeding the given value.
     */
    public void setSize(int n) {
    	if(n < 0 || n >= Integer.MAX_VALUE) throw new IllegalArgumentException("Inappropriate size specified for lattice: "+n);
		delayUpdate = true;
		int i = 0;
		//set lattice dimensions to value that gives at least n sites  
    	for(i=1; i<=n; i++) if(factory.crystal.getBasis().size() * factory.space.powerD(i) >= n) break;
    	setDimensions(factory.space.makeArrayD(i));//does a rebuild() call, but is blocked by delayUpdate
		((Factory)creator()).setDimensions(factory.space.makeArrayD(i));//does a rebuild() call, but is blocked by delayUpdate
		delayUpdate = false;
		rebuild();//also does update()
    }
	/**
	 * Performs same actions as setSize(int), then size of primitive vectors are
	 * adjusted such that lattice will fit in the given dimensions.  Assumes all
	 * dimension values are equal (future development will address case of non-
	 * square dimensions).
	 */
	public void setSize(int n, double[] dimensions) {
		if(n < 0 || n >= Integer.MAX_VALUE) throw new IllegalArgumentException("Inappropriate size specified for lattice: "+n);
		delayUpdate = true;
		int i = 0;
		//set lattice dimensions to value that gives at least n sites  
		for(i=1; i<=n; i++) if(factory.crystal.getBasis().size() * factory.space.powerD(i) >= n) break;
		setDimensions(factory.space.makeArrayD(i));//does a rebuild() call, but is blocked by delayUpdate
		((Factory)creator()).setDimensions(factory.space.makeArrayD(i));//does a rebuild() call, but is blocked by delayUpdate
    	//size primitive vectors to given dimensions
    	double[] size = new double[dimensions.length];
    	for(int j=0; j<dimensions.length; j++) size[j] = dimensions[j]/(double)i;
    	factory.crystal.getPrimitive().setSize(size); //does an update() call, but is blocked by delayUpdate
    	delayUpdate = false;
    	rebuild();//also does update()
    }
    
    /**
     * Translates the lattice so that the first site is positioned at the origin,
     */
    public void shiftFirstToOrigin() {
        Space.Vector shift = Space.makeVector(D);
        shift.ME(siteList.getFirst().coord.position());
        this.coord.translateBy(shift);
    }

    /** 
     * Sets the primitive and site factory to those specified by the given
     * Crystal class, and rebuilds the lattice.
     */
    public void setCrystal(Crystal crystal) {
        ((AtomTreeNodeGroup)node).removeAllChildren();//need to do this before setting new leaf factory, in case new factory has different groupFactory() field
        factory.crystal = crystal;
        ((Factory)creator()).setLeafFactory(crystal.getBasis().atomFactory());
        rebuild();
    }
    
    //do not want to allow primitive to be set independently of basis.  use setCrystal
//    /**
//     * Sets the primitive for this lattice to the one given, and 
//     * updates the site positions.
//     */
//    public void setPrimitive(Primitive primitive) {
//        this.primitive = primitive;
//        primitive.setLattice(this);
//        update();
//    }
    /**
     * Returns the primitive object used to construct this lattice.
     * Note that if the primitive is modified, changes will not be
     * reflected in this lattice until the update() method is called.
     */
    public Primitive getPrimitive() {return factory.crystal.primitive;}
    
    /**
     * Returns the atom type instance held by all atoms at the depth d.
     */
    protected AtomType getAtomType(int d) {
        if(d < 0 || d >= D) throw new IllegalArgumentException("Error in BravaisLattice.getFactory: inappropriate dimension specified");
        Atom a = this;
        int i = d;
        while(i < d) {a = ((AtomTreeNodeGroup)a.node).firstChildAtom(); i++;}
        return a.type;
    }
    
    //not carefully implemented
    public Atom nearestSite(Space.Vector r) {
    	AtomIteratorList iterator = new AtomIteratorList(siteList);
    	iterator.reset();
    	Atom site = null;
    	double r2min = Double.MAX_VALUE;
    	while(iterator.hasNext()) {
    		Atom a = iterator.nextAtom();
    		double r2 = a.coord.position().Mv1Squared(r);
    		if(r2 < r2min) {
    			r2min = r2;
    			site = a;
    		}
    	}
    	return site;
//        for(int i=D-1; i>=0; i--) {
//            idx[i] = (int)(r.x(i)/primitiveVectorLength[i] + 0.5);
//        }
//        return site(idx);
    }
    
    /**
     * Returns the site corresponding to the given index, applying periodic
     * boundaries if the index describes a site beyond the lattice.
     */
     //needs work to improve probable inefficiency with list.get
     //have not checked handling of values out of size of lattice
    public Atom site(int[] idx) {
        Atom site = this;
        for(int i=0; i<idx.length && !site.node.isLeaf(); i++) {
            int k = idx[i];
            while(k < 0) k += factory.dimensions[i];
            while(k >= factory.dimensions[i]) k -= factory.dimensions[i];
            site = ((AtomTreeNodeGroup)site.node).childList.get(k);
        }
      //  int i=0;
      //  do site = ((AtomTreeNodeGroup)site.node).childList.get(idx[i++]);
      //  while(!site.node.isLeaf() && i<idx.length);
        return site;
    }
    
    /**
     * Returns a list of all the Bravais lattice sites for this lattice.
     * Sites making up the lattice basis are not returned individually
     * if the basis is multiatomic.
     */
    public AtomList siteList() {return siteList;}
    
    /**
     * Number of sites in the lattice, equal to the size of siteList.
     * @return int number of sites
     */
    public int siteCount() {return siteList.size();}
    /**
     * Returns a list of all the leaf lattice sites for this lattice.  For a
     * one-atom basis, this is the same as the siteList; otherwise, returns
     * a list of all atoms at the lowest hierarchy in the basis, for all basis
     * sites.  Creates list on-the-fly from siteList, so will be expensive if
     * used repeatedly.
     * @return AtomList list of all leaf atoms
     */
    //developed for use by ConfigurationLattice
//    public AtomList leafList() {
//    	AtomIteratorTree leafIterator = new AtomIteratorTree();
//    	leafIterator.setLeafIterator();
//    	AtomList leafList = new AtomList();
//    	AtomIteratorList siteIterator = new AtomIteratorList(siteList);
//    	siteIterator.reset();
//    	while(siteIterator.hasNext()) {
//    		leafIterator.setBasis(siteIterator.next());
//    		while(leafIterator.hasNext()) leafList.add(leafIterator.next());
//    	}    	
//    	return leafList;
//    }
    /**
     * Returns the spatial dimension of this lattice.
     */
    public int D() {return D;}
    
    /**
     * Returns the event manager that registers listeners and notifies them
     * of events indicating changes in this lattice.
     * Part of AbstractLattice interface.
     */
    public SimulationEventManager eventManager() {return eventManager;}
    
    /**
     * Reconstructs all sites of the lattice.  Invoked when dimensions of
     * lattice are changed.
     */
    public void rebuild() {
    	if(delayUpdate) return;
        ((AtomTreeNodeGroup)node).removeAllChildren();
        creator().build(this);
        eventManager.fireEvent(rebuildEvent);
    }

    /**
     * Overrides the superclass method to return a simple string.
     * Terminates the chain of parentGroup calls that define the signature
     * of an atom in this lattice.
     */
    public String signature() {return "BravaisLattice";}
    
    /**
     * Causes all coordinates to update their position vectors, and
     * notifies any observers that a change has occurred.
     * Invoked by setPrimitive.  If primitive vectors are altered externally,
     * the modifying program must call this method to effect any changes
     * in the lattice coordinates.
     */
    public void update() {
    	if(delayUpdate) return;
        Space.Vector[] pVectors = factory.crystal.getPrimitive().vectors();
        Atom a = this;
        for(int i=0; i<D; i++) {
           ((ConfigurationLinear)a.creator().getConfiguration()).setOffset(pVectors[i]);
           primitiveVectorLength[i] = Math.sqrt(pVectors[i].squared());
           a = ((AtomTreeNodeGroup)a.node).firstChildAtom();
        }
        type.creator().getConfiguration().initializePositions(this);
        eventManager.fireEvent(allSiteEvent);
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
    /**
     * Sets up the neighbor list of each site according to the given criterion.
     */
    public void setupNeighbors(NeighborManager.Criterion criterion) {
        neighborCriterion = criterion;
        //inefficient method for setting up neighbors.  Repeats same calculation over all sites
//        AtomIteratorList iterator = new AtomIteratorList(siteList);
//        iterator.reset();
//        while(iterator.hasNext()) {
//            Site site = (Site)iterator.next();
//            site.neighborManager().setupNeighbors(siteList, criterion);
//        }
 
        //much more efficient approach:
        //set up neighbors for first site, then copy that neighbor structure to all other sites
        Site first = (Site)siteList.getFirst();
        int[] i0 = first.latticeCoordinate();
        first.neighborManager().setupNeighbors(siteList,criterion);
        
        //compute differences in index coordinates between each nbr and first site
        SiteIteratorNeighbor nbrIterator = new SiteIteratorNeighbor();
        nbrIterator.setBasis(first);
        nbrIterator.reset();
        int nbrCount = nbrIterator.size();
        int[][] delta = new int[nbrCount][D];
        int n = 0;
        int[] iN;
        while(nbrIterator.hasNext()) {            
            Site nbr = (Site)nbrIterator.next();
            if(nbr == first) {nbrCount--; continue;} //in case nbrIterator returns central atom itself
            iN = nbr.latticeCoordinate();
            for(int i=0; i<D; i++) delta[n][i] = iN[i] - i0[i];
            n++;
        }//end while
        
        //loop through all other atoms, adding neighbors identified by deltas
        //in index coordinates
        AtomIteratorList iterator = new AtomIteratorList(siteList);
        iterator.reset();
        iterator.next(); //skip first
        iN = new int[D];
        NbrCriterion nbrCriterion = new NbrCriterion();
        while(iterator.hasNext()) {
            nbrCriterion.root = (Site)iterator.nextAtom();
            i0 = nbrCriterion.root.latticeCoordinate();
            nbrCriterion.nbrList.clear();
            for(n=0; n<nbrCount; n++) {
                for(int i=0; i<D; i++) iN[i] = i0[i] + delta[n][i];
                nbrCriterion.nbrList.add(site(iN));
                //site.neighborManager().add to neighbor list
            }
            //could be made faster if we could add neighbors to neighborManager's list
            //in the appropriate order (pbc complicates doing this)
            nbrCriterion.root.neighborManager().setupNeighbors(siteList, nbrCriterion);
        }
        eventManager.fireEvent(resetNbrEvent);
    }//end of setupNeighbors
    private class NbrCriterion implements NeighborManager.Criterion {
        Site root;
        AtomList nbrList = new AtomList();
        public boolean areNeighbors(Site s0, Site s1) {
            return (s0 == root && nbrList.contains(s1))
                    || (s1 == root && nbrList.contains(s0));
        }
    }
    

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

 /**
  * Factory to construct an arbitrary-dimension Bravais lattice.
  */
  //need to update to permit modification of dimensions
  //make build(Atom) work using atom's parameters
protected static class Factory extends AtomFactoryTree {
    
     //dimension of lattice equals the number of primitive vectors
     //dimension of embedding space equals the length of each primitive vector
     //(trying to leave open possibility of having, e.g., a 2-D lattice in a 3-D space)
    Factory(Space space, Crystal crystal, int[] dimensions) {
		super(space, AtomSequencerSimple.FACTORY, crystal.getBasis().atomFactory(), dimensions, configArray(space, crystal.getPrimitive().vectors()));
		if(dimensions.length > space.D()) throw new IllegalArgumentException("Error: dimensions inconsistent with space");
		this.crystal = crystal;
		this.dimensions = (int[])dimensions.clone();
	}
//	//01/20/04 (DAK) deleted this constructor and replaced with the one above
//	Factory(Space space, AtomSequencer.Factory seqFactory, AtomFactory siteFactory, int[] dimensions, Primitive primitive) {
//		this(space, seqFactory, new Crystal(primitive, new Basis(space, 1, siteFactory)), dimensions);
//	}
    
    public Atom build(AtomTreeNodeGroup parent) {
        BravaisLattice group = new BravaisLattice(this.space, groupType, parent);
        return build(group);
    }
    
    public Atom build(Atom atom) {
        if(!(atom instanceof BravaisLattice)) throw new IllegalArgumentException("Error in BravaisLattice.Factory.build(Atom): Attempt to rebuild lattice from atom that is not a BravaisLattice instance");
        super.build(atom);
        BravaisLattice group = (BravaisLattice)atom;
        group.factory = this;
        AtomIteratorTree leafIterator = new AtomIteratorTree();
        leafIterator.setRoot(group);
        leafIterator.reset();
        group.siteList.clear();
        group.siteList.addAll(leafIterator);
        group.update();
        return group;
    }

    private static Configuration[] configArray(Space space, Space.Vector[] pVectors) {
        if(pVectors.length != space.D()) throw new IllegalArgumentException("Error in BravaisLattice.Factory constructor:  number of primitive vectors inconsistent with dimension of space");
        Configuration[] array = new Configuration[pVectors.length];
        for(int i=0; i<array.length; i++) {
            array[i] = new ConfigurationLinear(space);
            ((ConfigurationLinear)array[i]).setOffset(pVectors[i]);
        }
        return array;
    }
    
    void setDimensions(int[] dimensions) {
        if(dimensions.length > space.D()) throw new IllegalArgumentException("Incorrect size of array argument in BravaisLattice.Factory.setDimensions");
        System.arraycopy(dimensions, 0, this.dimensions, 0, dimensions.length);
        rootFactory().setNAtoms(this.dimensions);
    }
    int[] getDimensions() {return dimensions;}
    
    Crystal crystal;
    int[] dimensions;
}//end of Factory

    /**
     * Main method to demonstrate use of BravaisLattice and to aid debugging
     */
    public static void main(String[] args) {
        System.out.println("main method for BravaisLattice");
        Simulation sim = Simulation.instance;
        Space space = new Space2D();
        int D = space.D();
        PrimitiveOrthorhombic primitive = new PrimitiveOrthorhombic(space);
        final int nx = 4;
        final int ny = 5;
        BravaisLattice lattice = BravaisLattice.makeLattice(space, 
                                new Site.Factory(space),
                                new int[] {nx,ny},
                                primitive);        
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
        AtomActionAdapter printSites = new AtomActionAdapter() {
            public void actionPerformed(Atom s) {
                System.out.print(s.coord.position().toString()+" ");
       //         System.out.println(((Site)s).latticeCoordinate()[1]);
                if(((Site)s).latticeCoordinate()[1]==ny-1) System.out.println();
            }
        };
 //       iterator.allAtoms(printSites);
        System.out.println();
        
        System.out.println();
        System.out.println("Changing primitive vector");
        primitive.setSize(new double[] {1.0, 0.5});
        lattice.update();
 //       lattice.setPrimitiveVector(new Space.Vector[] {Space.makeVector(new double[] {0.,1.}),
 //                                                      Space.makeVector(new double[] {0.5,0.})});
 //       iterator.allAtoms(printSites);
        System.out.println();
        
        Atom testSite = lattice.site(new int[] {1,1});
        int[] idx = null;
        
        System.out.println();
        System.out.println("Translating lattice by (-1.0, 2.0)");
        lattice.coord.translateBy(Space.makeVector(new double[] {-1.0, 2.0}));
 //       iterator.allAtoms(printSites);
        System.out.println();
 
        System.out.println();
        System.out.println("Translating to origin");
        lattice.shiftFirstToOrigin();
 //       iterator.allAtoms(printSites);
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
 */       
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
//        nbrIterator.reset(IteratorDirective.BOTH);
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
 //       nbrIterator.reset(IteratorDirective.BOTH);
        while(nbrIterator.hasNext()) {  //print out coordinates of each site
            System.out.print(nbrIterator.next().toString()+" ");
        }
        System.out.println();
        System.out.println();
 /* */      
    }//end of main
}//end of BravaisLattice