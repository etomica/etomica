package etomica;

/**
 * Loops through molecule pairs given between two groups, formed
 * from 1 molecule and All its neighbors.  Molecule is specified itself
 * or in terms of one of its descendant atoms, via an iterator directive.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 08/04/02 Changed reset(atom) to use isDescendedFrom rather than just checking parent; also changed to check against iteratorDirective.direction
  *          Changes conducted as part of effort to resolve problems with speciesPistonCylinder; iterating over groups in which one has a molecule layer and the other doesn't
  * 12/06/02 (DAK) made class not final so could subclass in sufactant module 
  * 01/02/03 (DAK) reset checks if basis is null
  */
 
public class ApiIntergroup1A implements AtomPairIterator {
    
    protected Atom group1; 
    protected Atom group2;
    
    protected /*final*/ AtomIterator atomIterator;
    
    protected Atom referenceAtom;
    protected final IteratorDirective localDirective = new IteratorDirective(IteratorDirective.BOTH);
    protected final AtomPair pair;
    
    public ApiIntergroup1A(Simulation sim) {
        pair = new AtomPair(sim.space);
        atomIterator = sim.iteratorFactory.makeIntergroupNbrIterator();
    }
    
	public void all(AtomSet basis, IteratorDirective id, final AtomSetAction action) {
		 if(!(basis instanceof AtomPair && action instanceof AtomPairAction)) return;
		 all((AtomPair)basis, id, (AtomPairAction)action);
	}
	public void all(AtomPair basis, IteratorDirective id, AtomPairAction action) {
		Atom atom = id.atom1();
		if(atom == null || basis == null || action == null) return;
		Atom group1 = basis.atom1();//assume group1 preceeds group2
		Atom group2 = basis.atom2();
		if(group1 == group2) throw new IllegalArgumentException("Improper basis given to ApiIntergroup1A: Basis atoms must be different");

		AtomTreeNode referenceNode = atom.node.childWhereDescendedFrom(group1.node);
		if(referenceNode != null && localDirective.direction().doUp()) {
			action.wrapper.pair.atom1 = referenceNode.atom;
			atomIterator.all(group2, id.set(referenceNode.atom), action.wrapper);
		} else {
			referenceNode = atom.node.childWhereDescendedFrom(group2.node);
			if(referenceNode != null && localDirective.direction().doDown()) {
				action.wrapper.pair.atom1 = referenceNode.atom;
				atomIterator.all(group1, id.set(referenceNode.atom), action.wrapper);
			}
		}
	}//end all

   /**
     * Identifies the groups that are the parents of the atoms
     * to be given by the iterator.  The arguments must be different instances.
     * It is expected, but not verified, that a1.preceeds(a2).
     */
    public void setBasis(Atom a1, Atom a2) {
        if(a1 == a2)
            throw new IllegalArgumentException("Improper basis given to ApiIntergroup1A");
        group1 = a1;
        group2 = a2;
    }
    
    /**
     * Returns the number of pairs capable of being given by this iterator, based
     * on the most recent specification of the atom given to reset.
     */
    public int size() {return atomIterator.size();}        
    
    public boolean hasNext() {return atomIterator.hasNext();}
    
    public void reset(IteratorDirective id) {
        localDirective.set(id.direction());
        reset(id.atom1());
    }
    
    /**
     * Resets the iterator using the current reference atom, or unsets if no
     * atom was previously given.
     */
    public void reset() {
        if(referenceAtom == null || group1 == null || group2 == null) atomIterator.unset();
        else atomIterator.reset(localDirective.set(IteratorDirective.BOTH).set(referenceAtom));
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom.
     */
    //may need rewrite to use childWhereDescendedFrom
    public void reset(Atom atom) {
        if(atom == null || group1 == null || group2 == null) return;
        referenceAtom = atom;
        if(referenceAtom.node.isDescendedFrom(group1) && localDirective.direction().doUp()) {
            atomIterator.setBasis(group2);
        } else if(referenceAtom.node.isDescendedFrom(group2) && localDirective.direction().doDown()) {
            atomIterator.setBasis(group1);
        } else {
            atomIterator.unset();
            return;
        }
/*        Atom referenceGroup = referenceAtom.node.parentGroup();
        if(referenceGroup == group1) {
            atomIterator.setBasis(group2);
        } else if(referenceGroup == group2) {
            atomIterator.setBasis(group1);
        } else {
            atomIterator.setBasis(null);
            atomIterator.reset();
            return;
        }*/
        atomIterator.reset(localDirective.set(referenceAtom));
        pair.atom1 = referenceAtom;
    }
    
    public AtomPair next() {
        pair.atom2 = atomIterator.next();
        pair.reset();
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allPairs(AtomPairAction act) {
        if(referenceAtom == null) return;
        atomIterator.allAtoms(act.wrapper);
    }
    
}  //end of class AtomPairIterator
    
