package etomica;

/**
 * General iterator in which order of atoms is that given by the 
 * linked list of atoms, which changes only if atoms are added or removed to/from the phase.
 * All atoms are considered neighbors of each other.  Subclasses can vary in the range
 * of the atoms iterated, as specified by the defaultFirstAtom and defaultLastAtom methods.
 * Phase.AtomIterator subclasses this and defines first/last atoms of iterator to coincide
 * with first/last atoms of atom sequence in the phase.  Species.AtomIterator subclasses
 * to include only the atoms in the species, and so on.
 *
 * @author David Kofke
 */
public class AtomIteratorSequential extends AtomIteratorAbstract  {

    /**
     * Initial state is hasNext = false.
     */
    public AtomIteratorSequential() {
        this(false);
    }
    public AtomIteratorSequential(boolean isLeafIterator) {
        super(isLeafIterator);
    }
    public AtomIteratorSequential(Atom a) {
        this(a, false);
    }
    //instead make an enumerated type to key for leaf iterator
    public AtomIteratorSequential(Atom a, boolean isLeafIterator) {
        super(isLeafIterator);
        setBasis(a);
    }
    
    
    public Atom firstUpNeighbor(Atom a) {
        Atom first = defaultFirstAtom();
        if(a.preceeds(first)) return first;
        else if(a.preceeds(defaultLastAtom())) return a.nextAtom();
        else return null;
    }
    
    public Atom firstDownNeighbor(Atom a) {
        Atom last = defaultLastAtom();
        if(last == null || last.preceeds(a)) return last;
        else if(defaultFirstAtom().preceeds(a)) return a.previousAtom();
        else return null;
    }
    
    public boolean isOrdered(Atom atom1, Atom atom2) {
        if(atom1 == null || atom2 == null || !contains(atom1) || !contains(atom2)) return false;
        else return atom1 == atom2 || atom1.preceeds(atom2);
    }

    public Atom next() {
        Atom nextAtom = atom;
        //for debugging
    /*    if(atom == null) {
            System.out.println("error:  unexpected null in atomiteratorsequential");
        }*/
        //------------
        atom = upListNow ? atom.nextAtom() : atom.previousAtom();
        hasNext = nextAtom != terminator;
        if(!hasNext && upListNow && doGoDown) {//done going up and now prepare to go down
            atom = isNeighborIterator ? firstDownNeighbor(setAtom) : setAtom.previousAtom();
            terminator = terminator2;
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
        for(Atom a=defaultFirstAtom(); a!=null; a=a.nextAtom()) {act.actionPerformed(a);}
    }
    
    //main method to demonstrate and test this class
    public static void main(String[] args) throws java.io.IOException {
        
        java.io.BufferedReader in = new java.io.BufferedReader(new java.io.InputStreamReader(System.in));
        Simulation.instance = new Simulation();
        Phase phase = new Phase();
        Phase phase2 = new Phase();
        Species species = new SpeciesSpheres();
        species.setNMolecules(8);
        Simulation.instance.elementCoordinator.go();
        AtomIterator iterator = phase.iteratorFactory().makeAtomIterator();
        iterator.reset();
        Atom atom3 = null;
        Atom atom6 = null;
        int i=0;
        while(iterator.hasNext()) {
            Atom next = iterator.next();
            System.out.println(next.toString());
            if(i == 3) atom3 = next;
            if(i == 6) atom6 = next;
            i++;
        }
        IteratorDirective directive = new IteratorDirective();
        
        String line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(3)UP: 3..7/4..7");
        iterator.reset(directive.set(atom3).set(IteratorDirective.UP));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());
        
        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("()DOWN; 7..0/6..0");
        iterator.reset(directive.set().set(IteratorDirective.DOWN));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(3)DOWN; 3..0/2..0");
        iterator.reset(directive.set(atom3).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(3)NEITHER; 3/null");
        iterator.reset(directive.set(atom3).set(IteratorDirective.NEITHER));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(3)BOTH; 3..7,2..0/4..7,2..0");
        iterator.reset(directive.set(atom3).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(0)UP; 0..7/1..7");
        iterator.reset(directive.set(phase.firstAtom()).set(IteratorDirective.UP));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(0)DOWN; 0/null");
        iterator.reset(directive.set(phase.firstAtom()).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(0)BOTH; 0..7/1..7");
        iterator.reset(directive.set(phase.firstAtom()).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(7)BOTH; 7..0/6..0");
        iterator.reset(directive.set(phase.lastAtom()).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(X)BOTH; null/null");
        iterator.reset(directive.set(phase2.firstAtom()).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(X)UP; null/null");
        iterator.reset(directive.set(phase2.firstAtom()).set(IteratorDirective.UP));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        Species species2 = new SpeciesSpheres();
        species.setNMolecules(4);
        species2.setNMolecules(5);
        Simulation.instance.elementCoordinator.go();
        
        System.out.println("second species added");
        
        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(0)UP; 0..N/1..N");
        iterator.reset(directive.set(phase.firstAtom()).set(IteratorDirective.UP));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(0)DOWN; 0/null");
        iterator.reset(directive.set(phase.firstAtom()).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());
 
        iterator = phase.makeAtomIterator(species2);
        System.out.println("iterator changed to second-species iterator");
        
        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(0)UP; null/0..N,1");
        iterator.reset(directive.set(phase.firstAtom()).set(IteratorDirective.UP));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(0)DOWN; null/null");
        iterator.reset(directive.set(phase.firstAtom()).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine(); iterator.setAsNeighbor(false);
        System.out.println("(N)DOWN; N..0/N-1..0,1");
        iterator.reset(directive.set(phase.lastAtom()).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) System.out.println(iterator.next().toString()); iterator.setAsNeighbor(true);
        iterator.reset(directive); System.out.println(); while(iterator.hasNext()) System.out.println(iterator.next().toString());

        line = in.readLine();  System.exit(0);
    }//end of main
        
} //end of AtomIteratorSequence
