package etomica;

/**
 * Builds an atom group that comprises atoms arranged in a tree structure.  The
 * atom given by the build() method will be the root of the tree.  It will have
 * a pre-set number of child atoms that are leaf atoms of the tree,
 * or are themselves the root of a subtree.  The leaf atoms that terminate
 * this hierarchy are built by the specified atom factory.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 09/04/02 (DAK) added leafFactory method
  * 09/16/02 (DAK) added parentFactory, rootFactory methods; modified constructors so
  *                that parent factory is set when child AtomFactoryTree is made
  * 09/26/02 (DAK) added line to setLeafFactory to update group.childrenAreGroups
  */
  
public class AtomFactoryTree extends AtomFactoryHomo {
    
    /**
     * Constructs atom tree.
     *
     * @param space Space that makes an instance of a coordinate for each atom.
     * @param leafFactory  Atom factory that builds the leaf atoms of the tree.
     * @param nAtoms Array specifiying the number of atoms at each level in the tree.
     *               nAtoms[0] gives the number of child atoms for the root, nAtoms[1] the number of
     *               child atoms under each of the root's child atoms, etc.
     */
    public AtomFactoryTree(Simulation sim, AtomFactory leafFactory, int[] nAtoms) {
        super(sim, subFactory(sim, leafFactory, nAtoms, null), nAtoms[0]);
        if(childFactory() != leafFactory) ((AtomFactoryTree)childFactory()).parentFactory = this;
        depth = nAtoms.length;
    }
    
    /**
     * Constructor that provides an array to specify the configuration
     * for each level in the tree.  config[0] arranges the child atoms of the
     * root, config[1] the child atoms of each root's child, etc.  Note that
     * arranging the atoms at each level implies arranging the center of mass
     * of all the atoms below it in the tree, and thus prescribes in part the positioning of 
     * all atoms below its level.
     */
    public AtomFactoryTree(Simulation sim, AtomFactory leafFactory, int[] nAtoms, 
                            Configuration[] config) {
        super(sim, subFactory(sim, leafFactory, nAtoms, config), 
                nAtoms[0], BondInitializer.NULL, config[0]);
        if(childFactory() != leafFactory) ((AtomFactoryTree)childFactory()).parentFactory = this;
        depth = nAtoms.length;
    }
    
    //method used by constructor to determine the child factory
    private static AtomFactory subFactory(Simulation sim, AtomFactory leafFactory, 
                                          int[] nAtoms, Configuration[] config) {
        if(nAtoms.length < 1) throw new IllegalArgumentException("Error: Attempt to prescribe zero-dimensional lattice in AtomFactoryLattice" );
        if(config != null && nAtoms.length != config.length) throw new IllegalArgumentException("Error: incompatible specification of nAtoms and config in AtomFactoryTree constructor");
        if(nAtoms.length == 1) return leafFactory; 
        else {
            int[] newDim = new int[nAtoms.length-1];//arraycopy
            for(int i=1; i<nAtoms.length; i++) newDim[i-1] = nAtoms[i];
            if(config == null) {
                return new AtomFactoryTree(sim, leafFactory, newDim);
            } else {
                Configuration[] newConfig = new Configuration[config.length-1];//arraycopy
                for(int i=1; i<config.length; i++) newConfig[i-1] = config[i];
                return new AtomFactoryTree(sim, leafFactory, newDim, newConfig);
            }
        }
    }//end of subFactory
    
    public void setNAtoms(int[] n) {
        if(n.length != depth) throw new IllegalArgumentException("AtomFactoryTree.setNAtoms(int[]) attempt to set value inconsistent with depth specified when factory was instantiated");
        AtomFactoryHomo factory = this;
        for(int i=0; i<n.length; i++) {
            factory.setAtomsPerGroup(n[i]);
            if(i < n.length -1) factory = (AtomFactoryHomo)factory.childFactory();
        }
    }
    
/*    private void setParentFactory(AtomFactoryTree parent) {
        parentFactory = parent;
    }*/
    public AtomFactoryTree parentFactory() {return parentFactory;}
    public AtomFactoryTree rootFactory() {
        return (parentFactory == null) ? this : parentFactory.rootFactory();
    }
    
    /**
     * Returns the factory that produces the leaf atoms of the tree.
     */
     public AtomFactory getLeafFactory() {
        return (childFactory() instanceof AtomFactoryTree) ?
                ((AtomFactoryTree)childFactory()).getLeafFactory() : childFactory();
     } 
     
     /**
      * Sets the factory that makes the leaf atoms of the tree.
      */
     public void setLeafFactory(AtomFactory factory) {
        if(childFactory instanceof AtomFactoryTree) ((AtomFactoryTree)childFactory).setLeafFactory(factory);
        else {
            childFactory = factory;
            groupType.childrenAreGroups = factory.isGroupFactory();// 09/26/02 (DAK)
        }     
     }
    
    //number of layers of atoms below the root atom
    int depth;
    private AtomFactoryTree parentFactory;
    
    public static void main(String[] args) {
        Simulation sim = new Simulation();
        AtomFactoryMono leafFactory = new AtomFactoryMono(sim);
        int[] nA = new int[] {2, 1, 3};
        AtomFactoryTree treeFactory = new AtomFactoryTree(sim, leafFactory, nA);
        Phase phase = new Phase();
        Species species = new Species(sim, treeFactory);
        species.setNMolecules(1);
        Atom atom = treeFactory.makeAtom();
        sim.elementCoordinator.go();
        
        AtomIteratorTree iterator = new AtomIteratorTree();
        iterator.setBasis(phase.speciesMaster);
        for(int i=0; i<6; i++) {
            System.out.println("i = "+i);
            iterator.setIterationDepth(i);
            iterator.reset();
            while(iterator.hasNext()) System.out.print(iterator.next().toString());
            System.out.println();
        }
        
        treeFactory.setNAtoms(new int[] {2,2,1});
        AtomIteratorMolecule moleculeIterator = new AtomIteratorMolecule(phase);
        while(moleculeIterator.hasNext()) treeFactory.build(moleculeIterator.next());
        for(int i=0; i<6; i++) {
            System.out.println("i = "+i);
            iterator.setIterationDepth(i);
            iterator.reset();
            while(iterator.hasNext()) System.out.print(iterator.next().toString());
            System.out.println();
        }
    }
    
}//end of AtomFactoryTree
