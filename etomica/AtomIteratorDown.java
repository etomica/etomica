package etomica;        
/**
 * Iterator that progresses down a list of atoms.
 * Order of atoms is that given by the linked list of atoms, which changes
 * only if atoms are added or removed from the phase.
 *
 * @author David Kofke
 */
public class AtomIteratorDown extends AtomIteratorAbstract  {
    
    private final Phase phase;

    /**
     * Sets phase but does not reset iterator.  
     * Initial state is hasNext = false.
     * Default initiation flag is INCLUDE_FIRST.
     */
    public AtomIteratorDown(Phase p) {
        this(p, AtomIterator.INCLUDE_FIRST);
    }
    public AtomIteratorDown(Phase p, AtomIterator.Initiation init) {
        super(init);
        phase = p; 
    }

    /**
     * Phase in which this iterator operates.
     */
    public Phase phase() {return phase;}

    /**
     * Returns the last atom in the phase, which is the first atom
     * returned by this iterator.
     */
    public Atom defaultFirstAtom() {return phase.lastAtom();}
    
    /**
     * Returns the first atom in the phase, which is the default
     * last atom returned by this iterator.
     */
    public Atom defaultLastAtom() {return phase.firstAtom();}

    /**
     * Returns true if atom2 is downlist (preceeds) atom1.
     */
    public boolean isOrdered(Atom atom1, Atom atom2) {return atom2.preceeds(atom1);}
        
    public Atom next() {
        Atom nextAtom = atom;
        atom = atom.previousAtom();
        hasNext = nextAtom != terminator;
        return nextAtom;
    }
    
    /**
     * Performs the given action on all atoms in the phase, starting from the last to the first.
     */
    public void allAtoms(final AtomAction act) {
        for(Atom a=phase.lastAtom(); a!=null; a=a.previousAtom()) {act.actionPerformed(a);}
    }
} //end of AtomIteratorDown
