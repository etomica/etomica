package etomica;

/**
 * Atom iterator that traverses the elements of an AtomList.
 * Configurable to permit iteration up and/or down list, beginning
 * with any specified atom in list.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 09/01/02 (DAK) modified nextLinker method to properly handle case of NEITHER direction
  */
public final class AtomIteratorList extends AtomIterator {
    
    private AtomList list;
	private AtomTreeNodeGroup basis;
    
	private AtomLinker first;
	private AtomLinker.Tab header, terminator;
	private IteratorDirective.Direction direction = IteratorDirective.UP;
	private boolean skipFirstAtom = false;
	
    private boolean upListNow, doGoDown;
	private AtomLinker next;
    
    /**
     * Constructs a new iterator using an empty list as its basis for iteration.
     */
	public AtomIteratorList() {
	    this(AtomList.NULL);
	}
	/**
	 * Loops through the given iterator as currently reset and constructs a
	 * new list of atoms from its iterates, and constructs a new iterator
	 * using this new list as its basis for iteration.  Iterator is reset
	 * for iteration upon construction (constructor calls reset()).
	 */
	public AtomIteratorList(AtomIterator iterator) {
	    this(new AtomList(iterator));
	}
	/**
	 * Constructs a new iterator using the given list as its basis for iteration.
	 * Iterator is reset for iteration upon constrcution (constructor calls reset()).
	 */
	public AtomIteratorList(AtomList list) {
	    setBasis(list);
	    reset();
	}
	
	/**
	 * Makes a linked-list copy of the given list and constructs
	 * an iterator using the copy as a basis.  Useful to iterate over
	 * a list of atoms while doing operations that might change
	 * their order in the original list.
	 */
	public static AtomIteratorList makeCopylistIterator(AtomList list) {
	    return new AtomIteratorList(new AtomIteratorList(list));
	}
	
	
	/**
	 * Returns true if this iterator has another iterate, false otherwise.
	 */
	public boolean hasNext() {return next.atom != null;}
	
	public void setSkipFirstAtom(boolean b) {skipFirstAtom = b;}
	public boolean isSkipFirstAtom() {return skipFirstAtom;}
	
	/**
	 * Sets the childList of the given atom as the basis for iteration.
	 * If atom is a leaf, sets empty list as basis.  Sets given atom to be
	 * that returned by the getBasis method.
	 */
    public void setBasis(Atom atom){
        if(atom == null || atom.node.isLeaf()) {
            setBasis(AtomList.NULL);
            return;
        }
        setBasis((AtomTreeNodeGroup)atom.node);
    }
    
	/**
	 * Sets the childList of the given node as the basis for iteration.
	 * Sets node's atom to be that returned by the getBasis method.
	 */
    public void setBasis(AtomTreeNodeGroup node) {
        basis = node;
        this.list = basis.childList;
        header = list.header;
        next = terminator = header;
    }   
    
    /**
     * Sets the given list of atoms as the basis for iteration.  The atoms
     * returned by this iterator will be those from the given list.  A subsequent
     * call to one of the reset methods is required before iterator is ready 
     * for iteration (until then, hasNext is false).  Sets basis to null.
     */
    public void setBasis(AtomList list) {
        basis = null;
        if(list == null) list = AtomList.NULL;
        this.list = list;
        header = list.header;
        next = terminator = header;
    }
    /**
     * Returns current basis as set via setBasis(Atom), or null if setBasis(AtomList)
     * was most recently invoked.
     */
    public Atom getBasis() {return (basis == null) ? null : basis.atom;}

    /**
     * Resets the iterator using the current values of first, terminator, and direction.  
     * Called by all the reset methods after they have set these fields.  Also called
     * by allAtoms to prepare for iteration.
     */
    private Atom doReset() {
        if(first == header || first == null) next = header;
        else {
            upListNow = direction.doUp();
            doGoDown = direction.doDown();
            next = first;
            if(next.atom == null || skipFirstAtom) nextLinker();//nextLinker keeps iterating until entry
        }                                                       //with atom is found or terminator is reached
        return next.atom;
    }

	public void all(Atom basis, IteratorDirective id, final AtomActive action) {
		if(basis == null || basis.node.isLeaf() || action == null) return;

		final AtomList list = ((AtomTreeNodeGroup)basis.node).childList;
		final AtomLinker.Tab header = list.header;
		AtomLinker.Tab terminator = header;
		AtomLinker first = null;
		AtomLinker next = null;
		final IteratorDirective.Direction direction = id.direction();
		boolean upListNow = direction.doUp();
		boolean doGoDown = direction.doDown();
		if(id.atomCount() == 0) {
			first = upListNow ? header.next : (doGoDown ? header.previous : header);
			terminator = header;
		} else {
			Atom atom = id.atom1();
			AtomTreeNode referenceNode = atom.node.childWhereDescendedFrom(basis.node);
			if(referenceNode == null) {//atom not descended from basis
				boolean preceed = atom.seq.preceeds(basis);
//				if(preceed && upListNow) first = header.next;
//				else if(!preceed && doGoDown) first = header.previous;
//				else return;
				//following statement is equivalent in effect to the above if-else block
				first = (preceed && upListNow) ? header.next : (!preceed && doGoDown) ? header.previous : (AtomLinker)header;
				terminator = header;
			} else {//atom descended from basis
				first = referenceNode.atom.seq;
				terminator = header;
			}
		} 
		if(first == header) return;
		next = first;
		if(upListNow) {
			if(id.skipFirst) next = next.next;
			if(terminator == null) {//end loop when first tag is encountered
				while(next.atom != null) {
					action.actionPerformed(next.atom);
					next = next.next;
				}//end while
			} else {//end loop when terminator is encountered
				while(next.atom != null || next != terminator) {//first part of "or" may be omitted, but in most cases it gives true (and thus precludes evaluation of part after "or") and it is perhaps faster to evaluate than comparison with terminator
					if(next.atom != null) action.actionPerformed(next.atom);
					next = next.next;
				}//end while
			}//end else
			if(!doGoDown) return;//for NEITHER case, handled at end of method
		}//end if(upListNow)
        
		if(doGoDown) {
			if(id.skipFirst || upListNow) {//skip first down iterate, either because iterator is set to do so, or it was already given in up iteration
				next = first.previous;
//				if(next.atom == null) {//need to advance to find first entry
//					if(next == header || terminator == null) return;
//					upListNow = false;//set so that nextLinker() proceeds in proper direction
//					nextLinker();//find first non-null entry
//					if(next.atom == null) return; //none found
//				}//end if
			}//end if
                
			if(terminator == null) {//end loop when first tag is encountered
				while(next.atom != null) {
					action.actionPerformed(next.atom);
					next = next.previous;
				}
			} else {//end loop when header is encountered
				while(next.atom != null || next != header) {//first part of "or" may be omitted, but in most cases it gives true (and thus precludes evaluation of part after "or") and it is perhaps faster to evaluate than comparison with header
					if(next.atom != null) action.actionPerformed(next.atom);
					next = next.previous;
				}
			}
			return;
		}//end if(doGoDown)
        
		//if reaching here uplistNow and doGoDown are both false at the outset;
		//then direction == NEITHER, and loop over only current atom if it is not null
		if(next.atom != null) action.actionPerformed(next.atom);		
	}
    /**
     * Performs action on all atoms as prescribed in most recent call to reset.
     * Set of atoms for this method is same as that which would be given
     * by a hasNext/next loop.
     */
    public void allAtoms(AtomAction action){
        doReset();
        if(upListNow) {
            if(terminator == null) {//end loop when first tag is encountered
                while(next.atom != null) {
                    action.actionPerformed(next.atom);
                    next = next.next;
                }//end while
            } else {//end loop when terminator is encountered
                while(next.atom != null || next != terminator) {//first part of "or" may be omitted, but in most cases it gives true (and thus precludes evaluation of part after "or") and it is perhaps faster to evaluate than comparison with terminator
                    if(next.atom != null) action.actionPerformed(next.atom);
                    next = next.next;
                }//end while
            }//end else
            if(!doGoDown) return;//for NEITHER case, handled at end of method
        }//end if(upListNow)
        
        if(doGoDown) {
            if(skipFirstAtom || upListNow) {//skip first down iterate, either because iterator is set to do so, or it was already given in up iteration
                next = first.previous;
                if(next.atom == null) {//need to advance to find first entry
                    if(next == header || terminator == null) return;
                    upListNow = false;//set so that nextLinker() proceeds in proper direction
                    nextLinker();//find first non-null entry
                    if(next.atom == null) return; //none found
                }//end if
            }//end if
                
            if(terminator == null) {//end loop when first tag is encountered
                while(next.atom != null) {
                    action.actionPerformed(next.atom);
                    next = next.previous;
                }
            } else {//end loop when header is encountered
                while(next.atom != null || next != header) {//first part of "or" may be omitted, but in most cases it gives true (and thus precludes evaluation of part after "or") and it is perhaps faster to evaluate than comparison with header
                    if(next.atom != null) action.actionPerformed(next.atom);
                    next = next.previous;
                }
            }
            return;
        }//end if(doGoDown)
        
        //if reaching here uplistNow and doGoDown are both false at the outset;
        //then direction == NEITHER, and loop over only current atom if it is not null
        if(next.atom != null) action.actionPerformed(next.atom);
        
    }//end of allAtoms
    
    /**
     * Sets iterator so that it is ready to go upList its entire list of iterates.
     */
    public Atom reset() {
        this.first = header.next;
        this.terminator = header;
        this.direction = IteratorDirective.UP;
        return doReset();
    }
    
    /**
     * Resets to begin with the i_th atom, iterating uplist to the end of list.
     * e.g., if i = 3, begins with the third atom in list.
     */
	public Atom reset(int i) {
	    this.first = list.entry(i);
	    this.terminator = header;
	    this.direction = IteratorDirective.UP;
	    return doReset();
	}

    /**
     * Resets using direction and atom specified in directive.  If directive indicates
     * no atom, iterates in given direction as described in reset(Direction) method.
     * If an atom is given by directive, iteration is as follows depending on direction
     * indicated by directive:<ul>
     *    <li>UP:      Proceeds up, starting with the given atom.
     *    <li>DOWN:    Proceeds down, starting with the given atom.
     *    <li>BOTH:    Proceeds up, starting with the given atom; then down, starting with the atom
     *                 preceding the given one.
     *    <li>NEITHER: Returns given atom and expires (singlet iteration)
     * </ul>
     */
     //should atom.seq be used for reset?
    public Atom reset(IteratorDirective id){
        switch(id.atomCount()) {
            case 0:  return reset(id.direction()); 
            case 1:  return reset(id.atom1(), id.direction()); 
            default: next = header; 
            return null;
        }
    }
    
    /**
     * Resets to iterate in the given direction.  If UP or BOTH, iterates from first to 
     * last in list; if DOWN, iterates from last to first in list; if NEITHER, no iteration 
     * is performed.
     */
    public Atom reset(IteratorDirective.Direction direction) {
        this.first = direction.doUp() ? header.next : (direction.doDown() ? header.previous : null);
        this.terminator = header;
        this.direction = direction;
        return doReset();
    }
    
    /**
     * Resets to begin with the given atom linker, proceeding upList to end.  
     * Does not check that the linker is an iterate of this iterator.
     */
    public Atom reset(AtomLinker first) {
        this.first = first;
        this.terminator = header;
        this.direction = IteratorDirective.UP;
        return doReset();
    }
    
    /**
     * Resets to begin with the given atom linker, proceeding in the given direction.  
     * Does not check that the linker is an iterate of this iterator.
     */
    public Atom reset(AtomLinker first, IteratorDirective.Direction direction) {
        this.first = first;
        this.terminator = header;
        this.direction =  direction;
        return doReset();
    }
    
    /**
     * Resets in reference to the given atom.  Finds the atom in the list and
     * calls reset(AtomLinker) in reference to its linker.  If atom is not in list,
     * hasNext will be false.  Iteration proceeds upList to the end.
     */
    public Atom reset(Atom atom) {
        this.first = list.entry(atom);
        this.terminator = header;
        this.direction = IteratorDirective.UP;
        return doReset();
    }

    /**
     * Resets in reference to the given atom.  Finds the atom in the list and
     * calls reset(AtomLinker) in reference to its linker.  If atom is not in list,
     * hasNext will be false.  Iterator proceeds in the given direction.
     */
    public Atom reset(Atom atom, IteratorDirective.Direction direction) {
        this.first = list.entry(atom);
        this.terminator = header;
        this.direction = direction;
        return doReset();
    }

    /**
     * Resets for new iteration, beginning with the atom of the first argument.
     * If first is an index, iterator is advanced to begin with the
     * next atom entry.
     */
    public Atom reset(AtomLinker first, AtomLinker.Tab terminator) {
        this.first = first;
        this.terminator = terminator;
        this.direction = IteratorDirective.UP;
        return doReset();
    }
    
    /**
     * Resets to begin iteration in the given direction, stopping when the specified tab
     * is encountered.  If direction is UP, iteration proceeds up list and ends when terminator
     * or header is encountered; likewise if direction is DOWN.  If direction is BOTH, 
     * proceeds up list from starting point until encountering header or terminator, 
     * and then down it from the starting point, again until encountering header or terminator.
     * If terminator is null, iteration halts in each direction when any tab (or the header) 
     * is encountered.<br>
     * To iterate completely in either or both directions (ignoring all tabs), use the 
     * reset(AtomLinker, Direction) method.<br>
     * Up iteration always begins by returning the given first linker
     * (unless it is a tab); down iteration always begins with the linker before the given
     * first one.
     */
    public Atom reset(AtomLinker first, AtomLinker.Tab terminator, IteratorDirective.Direction direction) {
        this.first = first;
        this.terminator = terminator;
        this.direction = direction;
        return doReset();
    }
    
    /**
     * Sets iterator such that hasNext() will return false.
     */
    public void unset() {
        this.first = header;
        this.direction = IteratorDirective.NEITHER;//unsets allAtoms method
        next = header;
        upListNow = doGoDown = false;
    }

    /**
     * Returns true if the given atom is in the list of iterates, false otherwise.
     */
	public boolean contains(Atom atom){
        return list.contains(atom);
	}
	
	/**
	 * Returns the total number of iterates that can be returned by this iterator, for
	 * its current list basis.
	 */
	public int size() {return list.size();}

	    
    public Atom next() {
        return nextLinker().atom;
    }
    
    /**
     * Returns the next atom in the list without advancing the iterator.
     */
    public Atom peek() {
        return next.atom;
    }
    
    public AtomLinker nextLinker() {
        AtomLinker nextLinker = next;
        next = upListNow ? next.next : (doGoDown ? next.previous : header);//9/1/02 added doGoDown? clause to handle NEITHER
        while(next.atom == null) {
            //if terminator is null we stop at the first encounter of a Tab linker
            //otherwise stop only if Tab linker is the specified terminator
            if(terminator == null || next == header || next == terminator) {
                if(upListNow && doGoDown) {//done going up and now prepare to go down
                    next = first.previous;
                    upListNow = false;
                } else {
                    break;
                }
            }
            else next = upListNow ? next.next : next.previous;
        }
        return nextLinker;
    }//end of nextLinker
    
    /**
     * Method to test and demonstrate use of class.
     */
    public static void main(String[] args) {
        
        Simulation sim = new Simulation();
        SpeciesSpheresMono species = new SpeciesSpheresMono();
        species.setNMolecules(10);//tested also for 0 and 1 molecule
        Phase phase = new Phase();
        sim.elementCoordinator.go();
//        AtomList atomList = phase.speciesMaster.atomList;
        AtomList atomList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
        
        AtomIteratorList iterator = new AtomIteratorList(atomList);
        Atom first = atomList.getFirst();
        Atom last = atomList.getLast();
        Atom middle = null;
        try {
            middle = atomList.get(atomList.size()/2);//exception thrown if list is empty
        } catch(IndexOutOfBoundsException ex) {}
        
        iterator.setSkipFirstAtom(true);
        IteratorDirective.testSuite(iterator, first, middle, last);
    }

}//end of AtomIteratorList

