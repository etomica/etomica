    package etomica;

/**
 * General iterator in which order of atoms is that given by the 
 * linked list of atoms, which changes
 * only if atoms are added or removed to/from the phase.
 *
 * @author David Kofke
 */
public abstract class AtomIteratorSequence extends AtomIteratorAbstract  {

    /**
     * Sets phase but does not reset iterator.  
     * Initial state is hasNext = false.
     */
    public AtomIteratorSequence() {
        super();
    }

    /**
     * Returns first atom in phase.
     */
//    public Atom defaultFirstAtom() {return phase.firstAtom();}
    
    /**
     * Returns last atom in phase.
     */
 //   public Atom defaultLastAtom() {return phase.lastAtom();}
    
    public Atom firstUpNeighbor(Atom a) {
        if(a.preceeds(defaultFirstAtom())) return defaultFirstAtom();
        else if(a.preceeds(defaultLastAtom())) return a.nextAtom();
        else return null;
    }
    
    public Atom firstDownNeighbor(Atom a) {
        if(defaultLastAtom().preceeds(a)) return defaultLastAtom();
        else if(defaultFirstAtom().preceeds(a)) return a.nextAtom();
        else return null;
    }
    
    public boolean isOrdered(Atom atom1, Atom atom2) {
        if(atom1 == null || atom2 == null || !contains(atom1) || !contains(atom2)) return false;
        else return atom1.preceeds(atom2);
    }

    public Atom next() {
        Atom nextAtom = atom;
        atom = upListNow ? atom.nextAtom() : atom.previousAtom();
        hasNext = nextAtom != terminator;
        if(!hasNext && upListNow && doGoDown) {//done going up and now prepare to go down
            atom = isNeighborIterator ? firstDownNeighbor(setAtom) : setAtom.previousAtom();
            terminator = terminator2;
            if(atom != null) atom = atom.previousAtom(); //uplist would have handled first atom, so skip it
            hasNext = (atom != null && nextAtom != terminator);
            upListNow = false;
        }
        return nextAtom;
    }

    /**
     * Performs the given action on all atoms in the phase.  Unaffected by any prior calls to reset.
     */
     
     //not implemented
    public void allAtoms(final AtomAction act) {
 //       for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
    }
} //end of AtomIteratorSequence
