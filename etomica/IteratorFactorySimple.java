package etomica;

public class IteratorFactorySimple implements IteratorFactory {
    
    public static final IteratorFactorySimple INSTANCE = new IteratorFactorySimple();
        
    public AtomIterator makeGroupIteratorSequential() {return new SequentialIterator(false);}//not neighbor iterator, so don't skip first atom

    public AtomSequencer makeSimpleSequencer(Atom atom) {return new AtomSequencerSimple(atom);}
        
    public AtomIterator makeIntragroupNbrIterator() {return new SequentialIterator(true);}//skip first atom for intragroup neighbor iteration
    public AtomIterator makeIntergroupNbrIterator() {return new SequentialIterator(false);}//don't skip for intergroup neighbor iteration
    
    public AtomSequencer makeNeighborSequencer(Atom atom) {return new AtomSequencerSimple(atom);}
    
    public Class simpleSequencerClass() {return AtomSequencerSimple.class;}
    
    public Class neighborSequencerClass() {return AtomSequencerSimple.class;}
    
    public AtomSequencer.Factory simpleSequencerFactory() {return AtomSequencerSimple.FACTORY;}
    
    public AtomSequencer.Factory neighborSequencerFactory() {return AtomSequencerSimple.FACTORY;}
   
/////////////////////////////////////////////////////////////////////////////////////////////
    /**
    * Iterates over all the children of a given atom group.
    *
    * @author David Kofke
    */

private static final class SequentialIterator implements AtomIterator {
    
    private AtomIteratorList listIterator = new AtomIteratorList();
    private AtomTreeNodeGroup basis;
    
    public SequentialIterator(boolean skipFirstAtom) {
        listIterator.setSkipFirstAtom(skipFirstAtom);
    }
    
    /**
     * Defines the atoms that are subject to iteration as the children of the
     * given atom.
     */
    public void setBasis(Atom a) {
        basis = (AtomTreeNodeGroup)a.node;
        listIterator.setBasis(basis);
    }
    
    /**
     * Resets iterator so that it will loop up the list of atoms beginning
     * from the first one.
     */
    public Atom reset() {return listIterator.reset();}
    
    /**
     * Sets to state in which hasNext is false.
     */
    public void unset() {listIterator.unset();}
        
    public boolean hasNext() {return listIterator.hasNext();}
    
    public boolean contains(Atom atom) {return listIterator.contains(atom);}
    
    /**
     * Resets for iteration according to the given directive.  If the directive does
     * not specify an atom, this is the same as the reset() method, except that the
     * direction of iteration is as given by the directive.  If an atom is specified,
     * iteration begins with it and proceeds up or down list from there.
     */
    public Atom reset(IteratorDirective id) {
        Atom refAtom = id.atom1();
        if(refAtom.node.parentNode() == basis) {
            return listIterator.reset(id.atom1().seq, id.direction());
        } else if(id.direction()==IteratorDirective.BOTH) {
            return listIterator.reset();
        } else {
            boolean refFirst = refAtom.seq.preceeds(basis.atom());
            if(refFirst && id.direction()==IteratorDirective.UP 
                || !refFirst && id.direction()==IteratorDirective.DOWN && !(refAtom==basis.atom())) {
                    return listIterator.reset();
            }
        }
        return null;
     //   return listIterator.reset(id);
    }
    
    /**
     * Returns the next atom in the iteration sequence.  Assumes that hasNext is
     * true; calling when hasNext is false can lead to unpredictable results, and
     * may or may not cause an error or exception.
     */
    public Atom next() {return listIterator.next();}
    
    public void allAtoms(AtomAction act) {
        listIterator.allAtoms(act);
    }
    
    public Atom getBasis() {
        return basis.atom;
    }
    
    public int size() {
        return listIterator.size();
    }
    
    /**
     * Method to test SequentialIterator
     */
    public static void main(String args[]) {
        Simulation sim = new Simulation(new Space2D());

        IteratorFactory iteratorFactory = new IteratorFactorySimple();
        sim.setIteratorFactory(iteratorFactory);
        int nAtoms = 6;       
	    SpeciesSpheresMono speciesSpheres = new SpeciesSpheresMono();
	    speciesSpheres.setNMolecules(nAtoms);
	    Potential potential = new P2HardSphere();
	    Phase phase = new Phase();
	    IntegratorHard integrator = new IntegratorHard();
        integrator.setTimeStep(0.01);
        Controller controller = new Controller();
		Simulation.instance.elementCoordinator.go();
		System.out.println("Starting MD");
		integrator.setMaxSteps(100);
		integrator.initialize();
		integrator.run();
		System.out.println("Done");
		
		AtomIterator iterator = iteratorFactory.makeGroupIteratorSequential();
	    iterator.setBasis(phase.getAgent(speciesSpheres));
		Atom first = null;
		Atom middle = null;
		Atom last = null;
		iterator.reset();
		int k = 0;
		while(iterator.hasNext()) {
		    Atom a = iterator.next();
		    if(k == 0) first = a;
		    else if(k == nAtoms/2) middle = a;
		    else if(k == nAtoms-1) last = a;
		    k++;
		}
	    
	    IteratorDirective.testSuite(iterator, first, middle, last);
	    
    }//end of SequentialIterator.main

    }//end of Iterator
/////////////////////////////////////////////////////////////////////////////////////////////

}