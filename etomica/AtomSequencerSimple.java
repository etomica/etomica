package etomica;

public final class AtomSequencerSimple implements AtomSequencer {
    
    private Atom next;
    private Atom previous;
    private final Atom atom;
    
    public AtomSequencerSimple(Atom a) {
        atom = a;
    }
    
    public Atom nextAtom() {return next;}
    public Atom previousAtom() {return previous;}
    
    /**
     * Sets atom following this one in linked list, 
     * and sets this to be that atom's previous atom in list.
     * 
     * @param atom the next atom in the list
     */
    public void setNextAtom(Atom a) {
        next = a;
        if(a != null) a.seq.setPreviousAtom(this.atom);
    }
    
    public void setPreviousAtom(Atom a) {
        previous = a;
    }

    /**
     * Sets this atom's previous atom to be null
     * 
     * @see setNextAtom
     */
    public void clearPreviousAtom() {previous = null;}

    
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
}