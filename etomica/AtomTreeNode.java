package etomica;

public interface AtomTreeNode {
    
    public AtomGroup parentGroup();
    public Atom parentMolecule();
    public SpeciesAgent parentSpeciesAgent();
    public Species parentSpecies();
    public Phase parentPhase();
    public Simulation parentSimulation();
    
    public boolean isDescendedFrom(Atom a);
    
    public Atom atom();
    public void setAtom(Atom a);
    
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
 
    
    
}