package etomica;

/**
 * Iterator that progresses up a list of atoms.
 * Order of atoms is that given by the linked list of atoms, which changes
 * only if atoms are added or removed to/from the phase.
 */
public class AtomIteratorUp extends AtomIteratorAbstract  {

    private final Phase phase;
    /**
     * Sets phase but does not reset iterator.  
     * Initial state is hasNext = false.
     * Default initiation flag is INCLUDE_FIRST.
     */
    public AtomIteratorUp(Phase p) {
        this(p, AtomIterator.INCLUDE_FIRST);
    }
    public AtomIteratorUp(Phase p, AtomIterator.Initiation init) {
        super(init);
        phase = p; 
    }

    public Phase phase() {return phase;}
    
    /**
     * Returns first atom in phase.
     */
    public Atom defaultFirstAtom() {return phase.firstAtom();}
    
    /**
     * Returns last atom in phase.
     */
    public Atom defaultLastAtom() {return phase.lastAtom();}
    
    public boolean isOrdered(Atom atom1, Atom atom2) {return atom1.preceeds(atom2);}

    public Atom next() {
        Atom nextAtom = atom;
        atom = atom.nextAtom();
        hasNext = nextAtom != terminator;
        return nextAtom;
    }
    /**
     * Performs the given action on all atoms in the phase.  Unaffected by any prior calls to reset.
     */
    public void allAtoms(final AtomAction act) {
        for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
    }
} //end of AtomIteratorUp

