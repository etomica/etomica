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
    
    private static AtomFactory subFactory(Space space, AtomFactory leafFactory, int[] nAtoms) {
        if(nAtoms.length < 1) throw new IllegalArgumentException("Error: Attempt to prescribe zero-dimensional lattice in AtomFactoryLattice" );
        if(nAtoms.length == 1) return leafFactory; 
        else {
            int[] newDim = new int[nAtoms.length-1];
            for(int i=1; i<nAtoms.length; i++) newDim[i-1] = nAtoms[i];
            return new AtomFactoryTree(space, leafFactory, newDim);
        }
    }
    
    public AtomFactoryTree(Space space, AtomFactory leafFactory, int[] nAtoms) {
        super(space, subFactory(space, leafFactory, nAtoms), nAtoms[0]);
    }
    
}
