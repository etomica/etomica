package etomica;

/* History of changes
 * 07/29/02 (DAK) In SequentialIterator, changed reset(IteratorDirective) to handle null atom1()
 * 08/13/02 (DAK) In SequentialIterator, changed to allow basis to be leaf atom.  This introduces
 *               a singlet iterator and currentIterator to handle the general case.
 * 12/05/02 (DAK) Revised reset(IteratorDirective) so that it properly handles a reference atom
 *                this is descended from the basis, but not is a direct child of the basis.
 */

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
    private AtomIteratorSinglet singletIterator = new AtomIteratorSinglet();
    private AtomIterator currentIterator = AtomIterator.NULL;
    private AtomTreeNodeGroup basis;
    
    public SequentialIterator(boolean skipFirstAtom) {
        listIterator.setSkipFirstAtom(skipFirstAtom);
    }
 
	public void all(AtomSet basis, IteratorDirective id, final AtomSetAction action) {
		 if(!(basis instanceof Atom && action instanceof AtomAction)) return;
		 all((Atom)basis, id, (AtomAction)action);
	}
    
	public void all(Atom basis, IteratorDirective id, final AtomAction action) {
		if(basis == null || action == null) return;
		if(basis.node.isLeaf()) {
			singletIterator.all(basis, id, action);
		} else {
			listIterator.all(basis, id, action);
		}
	}
   
    /**
     * Defines the atoms that are subject to iteration as the children of the
     * given atom.
     */
    public void setBasis(Atom a) {
        if(!(a.node instanceof AtomTreeNodeGroup)) {//leaf atom
            //System.out.println("uh oh");
            basis = null;
           singletIterator.setBasis(a);
           currentIterator = singletIterator;
        } else { //atom group
           basis = (AtomTreeNodeGroup)a.node;
           listIterator.setBasis(basis);
           currentIterator = listIterator;
        }
    }
    
    /**
     * Resets iterator so that it will loop up the list of atoms beginning
     * from the first one.
     */
    public Atom reset() {return currentIterator.reset();}
    
    /**
     * Sets to state in which hasNext is false.
     */
    public void unset() {currentIterator.unset();}
        
    public boolean hasNext() {return currentIterator.hasNext();}
    
    public boolean contains(Atom atom) {return currentIterator.contains(atom);}
    
    /**
     * Resets for iteration according to the given directive.  If the directive does
     * not specify an atom, this is the same as the reset() method, except that the
     * direction of iteration is as given by the directive.  If an atom is specified,
     * iteration begins with it and proceeds up or down list from there.
     */
    public Atom reset(IteratorDirective id) {
        Atom refAtom = id.atom1();
        if(refAtom == null) {
            //cast used as temporary fix
            return ((AtomIteratorList)currentIterator).reset(id.direction());
        }
        //before adding descendedFrom tests, did only this and reset using seq
//        if(refAtom.node.parentNode() == basis) {

        //this does same as current implementation, but is perhaps less efficient
  /*      if(refAtom.node.isDescendedFrom(basis)) {
            return (refAtom.node.parentNode() == basis) ? 
                ((AtomIteratorList)currentIterator).reset(id.atom1().seq, id.direction())
               : ((AtomIteratorList)currentIterator).reset(id.direction());*/
               
        //If reference atom is child of basis, start iteration with it
        if(refAtom.node.parentNode() == basis) {//cast used as temporary fix
            ((AtomIteratorList)currentIterator).reset(id.atom1().seq, id.direction());
        //If descended from basis, but not a child, set for full iteration of children of basis
        } else if(refAtom.node.parentNode().isDescendedFrom(basis)) {
            // *** should modify so that reset is done from child of basis leading to refatom
            ((AtomIteratorList)currentIterator).reset(id.direction());
        //Reference is not descended from basis, but iteration proceeds over all children regardless
        } else if(id.direction()==IteratorDirective.BOTH) {
            return currentIterator.reset();
        //None of the above holds, need to look carefully at direction from reference atom to
        //basis, and compare to iteration direction
        } else {
            boolean refFirst = refAtom.seq.preceeds(basis.atom());
            if(refFirst && id.direction()==IteratorDirective.UP 
                || !refFirst && id.direction()==IteratorDirective.DOWN && !(refAtom==basis.atom())) {
                    return currentIterator.reset();
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
    public Atom next() {return currentIterator.next();}
    
    public void allAtoms(AtomAction act) {
        currentIterator.allAtoms(act);
    }
    
    public Atom getBasis() {
        return (basis == null) ? null : basis.atom;
    }
    
    public int size() {
        return currentIterator.size();
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