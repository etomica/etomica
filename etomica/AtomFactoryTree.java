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
    }
    
    /**
     * Constructor that provides an array to specifies the configuration
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
    
}//end of AtomFactoryTree
