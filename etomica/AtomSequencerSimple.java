package etomica;

public final class AtomSequencerSimple extends AtomSequencer {
    
    public AtomSequencerSimple(Atom a) {
        super(a);
    }
    
    /**
     * Called when position of atom is changed by translateBy or translateTo
     * methods.  In this sequencer this method performs no action.
     */
    public void moveNotify() {}
    
    /**
     * Called when parent of atom is changed by the node.  
     * In this sequencer this method performs no action.
     */
    public void setParentNotify(AtomTreeNodeGroup newParent) {}
    
    /**
     * Returns true if this atom preceeds the given atom in the atom sequence.
     * Returns false if the given atom is this atom, or (of course) if the
     * given atom instead preceeds this one.
     */
    public boolean preceeds(Atom a) {
        //want to return false if atoms are the same atoms
        if(a == null) return true;
        if(atom.node.parentGroup() == a.node.parentGroup()) return atom.node.index() < a.node.index();//works also if both parentGroups are null
        int thisDepth = atom.node.depth();
        int atomDepth = a.node.depth();
        if(thisDepth == atomDepth) return atom.node.parentGroup().seq.preceeds(a.node.parentGroup());
        else if(thisDepth < atomDepth) return this.preceeds(a.node.parentGroup());
        else /*if(this.depth > atom.depth)*/ return atom.node.parentGroup().seq.preceeds(a);
    }
    
    public static final AtomSequencer.Factory FACTORY = new AtomSequencer.Factory() {
        public AtomSequencer makeSequencer(Atom atom) {
            return new AtomSequencerSimple(atom);
        }
    };
}