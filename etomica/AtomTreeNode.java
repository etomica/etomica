package etomica;

/**
 * Interface for a node in the atom tree.  Atoms are arranged in a tree structure 
 * to define collections of molecules, individual molecules, atom groups and atoms.  
 * All of these elements are represented by the Atom class, and the structure that
 * forms the larger elements from the smaller ones is maintained by the node class
 * associated with each Atom.  The node for each atom holds information about its
 * parent and its children (if any).  Leaf nodes have no children, but they are defined
 * to identify themselves as their only child (this design makes certain types of
 * pair loops easier to perform).<br>
 *
 * (The node for) a SpeciesMaster instance is the root of the atom tree.  Directly below
 * it are instances of SpeciesAgent class, which are constructed each Phase by each Species.
 * Below the SpeciesAgent are instances of Atoms that are logically the molecules of the
 * system.  The molecules may be leaf atoms, ending the tree, or they may contain child 
 * atoms (or other groups of atoms) used to form the molecule.
 *
 * @author David Kofke
 */

public interface AtomTreeNode {
    
    public AtomGroup parentGroup();
    public AtomTreeNodeGroup parentNode();
    public Atom parentMolecule();
    public SpeciesAgent parentSpeciesAgent();
    public Species parentSpecies();
    public Phase parentPhase();
    public Simulation parentSimulation();
    
    public boolean isDescendedFrom(Atom a);
    
    public Atom atom();
    
    public Atom firstChildAtom();
    public Atom lastChildAtom();
    public Atom firstLeafAtom();
    public Atom lastLeafAtom();
    public Atom getAtom(int n);
    public Atom randomAtom();
    public Atom[] childAtomArray();
    
    public void addAtom(Atom a);
    public void removeAtom(Atom a);
    public Atom removeAll();
    
    public void addAtomNotify(Atom a);
    public void removeAtomNotify(Atom a);
    
    public int leafAtomCount();
    public int childAtomCount();
    public boolean childrenAreGroups();
    public int depth();
    public void setDepth(int d);
    /**
     * Integer assigned to this atom by its parent molecule.
     * Assigned during construction of atom.
     */
    public int index();
    public void setIndex(int i);
    
    public boolean isLeaf();
    
    public void setParentGroup(AtomGroup group);
 
    public interface Factory {
        public AtomTreeNode makeNode(Atom atom);
    }
    
}