    package etomica;

/**
 * General iterator in which order of atoms is that given by the 
 * linked list of atoms, which changes
 * only if atoms are added or removed to/from the phase.
 *
 * @author David Kofke
 */
public abstract class AtomIteratorSequential extends AtomIteratorAbstract  {

    /**
     * Sets phase but does not reset iterator.  
     * Initial state is hasNext = false.
     */
    public AtomIteratorSequential() {
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
    
    //main method to demonstrate and test this class
    public static void main(String[] args) throws java.io.IOException {
        
        java.io.BufferedReader in = new java.io.BufferedReader(new java.io.InputStreamReader(System.in));
        Simulation.instance = new Simulation();
        Phase phase = new Phase();
        Species species = new SpeciesDisks();
        Simulation.instance.elementCoordinator.go();
        AtomIterator iterator = phase.iteratorFactory().makeAtomIterator();
        iterator.reset();
        Atom atom5 = null;
        Atom atom10 = null;
        while(iterator.hasNext()) {
            Atom next = iterator.next();
            System.out.println(next.debugIndex);
            if(next.debugIndex == 5) atom5 = next;
            if(next.debugIndex == 10) atom10 = next;
        }
        String line = in.readLine();
        if(line.equals("n")) iterator.setAsNeighbor(true);
        iterator.reset(new IteratorDirective().set(atom5));
        while(iterator.hasNext()) System.out.println(iterator.next().debugIndex);
        iterator.setAsNeighbor(false);
        
        line = in.readLine();
        if(line.equals("n")) iterator.setAsNeighbor(true);
        iterator.reset(new IteratorDirective().set().set(IteratorDirective.DOWN));
        while(iterator.hasNext()) System.out.println(iterator.next().debugIndex);
        
        line = in.readLine();
        if(line.equals("n")) iterator.setAsNeighbor(true);
        iterator.reset(new IteratorDirective().set(atom5).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) System.out.println(iterator.next().debugIndex);

        line = in.readLine();
        System.exit(0);
    }
        
} //end of AtomIteratorSequence
