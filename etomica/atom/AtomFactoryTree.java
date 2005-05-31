package etomica.atom;

import etomica.AtomFactory;
import etomica.AtomType;
import etomica.AtomTypeGroup;
import etomica.Conformation;
import etomica.Space;

/**
 * Builds an atom group that comprises atoms arranged in a tree structure.  The
 * atom given by the build() method will be the root of the tree.  It will have
 * a pre-set number of child atoms that are leaf atoms of the tree,
 * or are themselves the root of a subtree.  The leaf atoms that terminate
 * this hierarchy are built by the specified atom factory.
 *
 * @author David Kofke
 */
   
public class AtomFactoryTree extends AtomFactoryHomo {

    /**
     * Constructs atom tree.
     *
     * @param space Governing space
     * @param parentType the AtomType of the parent of the atoms made by this factory
     * @param nAtoms Array specifiying the number of atoms at each level in the tree.
     *               nAtoms[0] gives the number of child atoms for the root, nAtoms[1] the number of
     *               child atoms under each of the root's child atoms, etc.
     */
    public AtomFactoryTree(Space space, AtomSequencerFactory seqFactory, AtomTypeGroup parentType,
                            int[] nAtoms) {
        this(space, seqFactory, parentType, nAtoms, null);
    }
    
    /**
     * Constructor that provides an array to specify the configuration
     * for each level in the tree.  config[0] arranges the child atoms of the
     * root, config[1] the child atoms of each root's child, etc.  Note that
     * arranging the atoms at each level implies arranging the center of mass
     * of all the atoms below it in the tree, and thus prescribes in part the positioning of 
     * all atoms below its level.
     */
    public AtomFactoryTree(Space space, AtomSequencerFactory seqFactory, AtomTypeGroup parentType,
                                int[] nAtoms, Conformation[] config) {
        super(space, seqFactory, parentType, nAtoms[0], (config != null) ? config[0] : Conformation.NULL);
        setChildFactory(subFactory(space, seqFactory, nAtoms, config));
        if(childFactory() != null) ((AtomFactoryTree)childFactory()).parentFactory = this;
        depth = nAtoms.length;
    }
    
    //method used by constructor to determine the child factory
    private AtomFactory subFactory(Space space, AtomSequencerFactory seqFactory,
                            int[] nAtoms, Conformation[] config) {
        if(config != null && nAtoms.length != config.length) throw new IllegalArgumentException("Error: incompatible specification of nAtoms and config in AtomFactoryTree constructor");
        if(nAtoms.length == 1) return null; 
        int[] newDim = new int[nAtoms.length-1];//arraycopy
        for(int i=1; i<nAtoms.length; i++) newDim[i-1] = nAtoms[i];
        if(config == null) {
            return new AtomFactoryTree(space, seqFactory, (AtomTypeGroup)atomType, newDim);
        }
        Conformation[] newConfig = new Conformation[config.length-1];//arraycopy
        for(int i=1; i<config.length; i++) newConfig[i-1] = config[i];
        return new AtomFactoryTree(space, seqFactory, (AtomTypeGroup)atomType, newDim, newConfig);
    }//end of subFactory
    
    public void setNAtoms(int[] n) {
        if(n.length != depth) throw new IllegalArgumentException("AtomFactoryTree.setNAtoms(int[]) attempt to set value inconsistent with depth specified when factory was instantiated");
        AtomFactoryHomo factory = this;
        for(int i=0; i<n.length; i++) {
            factory.setAtomsPerGroup(n[i]);
            if(i < n.length -1) factory = (AtomFactoryHomo)factory.childFactory();
        }
    }

    /**
     * Returns the factory above this one in the factory tree.
     */
    public AtomFactoryTree parentFactory() {return parentFactory;}
    
    /**
     * Returns the factory that is the root of the tree of factories, the one 
     * explicitly instantiated by the user (rather than one of the ones it instantiated
     * to make the tree). 
     */
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
        if (getLeafFactory() != null) throw new IllegalStateException("You can set the leaf factory only once!");
        if(childFactory instanceof AtomFactoryTree) {
            ((AtomFactoryTree)childFactory).setLeafFactory(factory);
        } else {
            childFactory = factory;
            childFactory.setSpecies(atomType.getSpecies());
        }     
     }
     
     /**
      * Returns the AtomType of the lowest-level non-null factory in the tree.
      * If a leaf factory has been assigned, this will be the type of that factory;
      * otherwise it will be the type of the AtomFactoryTree at the bottom of the
      * hierarchy.  This is needed to get the parent atom type when constructing the leaf factory.
      */
     public AtomType getLeafType() {
         if(childFactory == null) {
             return atomType;
         }
         if(childFactory instanceof AtomFactoryTree) {
             return ((AtomFactoryTree)childFactory).getLeafType();
         }
         return childFactory.getType();
     }
    
    //number of layers of atoms below the root atom
    int depth;
    private AtomFactoryTree parentFactory;
    
//    public static void main(String[] args) {
//        Simulation sim = new Simulation();
//        AtomFactoryMono leafFactory = new AtomFactoryMono(sim.space, AtomSequencerFactory.SIMPLE);
//        int[] nA = new int[] {2, 1, 3};
//        AtomFactoryTree treeFactory = new AtomFactoryTree(sim, leafFactory, nA);
//        Phase phase = new Phase(sim);
//        Species species = new Species(treeFactory);
//        species.setNMolecules(1);
//        phase.addMolecule(treeFactory.makeAtom(), species.getAgent(phase));
//        
//        AtomIteratorTree iterator = new AtomIteratorTree();
//        iterator.setRoot(phase.speciesMaster);
//        for(int i=0; i<6; i++) {
//            System.out.println("i = "+i);
//            iterator.setIterationDepth(i);
//            iterator.reset();
//            while(iterator.hasNext()) System.out.print(iterator.next().toString());
//            System.out.println();
//        }
//        
//        treeFactory.setNAtoms(new int[] {2,2,1});
//        AtomIteratorAllMolecules moleculeIterator = new AtomIteratorAllMolecules(phase);
//        while(moleculeIterator.hasNext()) treeFactory.build(moleculeIterator.nextAtom());
//        for(int i=0; i<6; i++) {
//            System.out.println("i = "+i);
//            iterator.setIterationDepth(i);
//            iterator.reset();
//            while(iterator.hasNext()) System.out.print(iterator.next().toString());
//            System.out.println();
//        }
//    }
    
}//end of AtomFactoryTree
