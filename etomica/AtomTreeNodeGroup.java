package etomica;

/**
 * Implementation of AtomTreeNode for non-leaf node.
 */

 /* History of changes
  * 08/21/02 (DAK) sketched out flatten method, but not implemented yet.
  * 08/12/03 (DAK) modified setParent method to notify children
  */
  
public class AtomTreeNodeGroup extends AtomTreeNode {
    
    public AtomTreeNodeGroup(Atom atom, AtomTreeNodeGroup parent) {
        super(atom, parent);
    }
        
    public int newChildIndex() {
        if(!childList.isEmpty()) { //siblings
            return (lastChildAtom().node.index()+1);
        }
        else {  //only child
            return(0);
        }
    }

	/**
	 * Invoke superclass method and then notify children.  Notification is
	 * needed, for example, to ensure their phase field is up to date.
	 * 
	 * @see etomica.AtomTreeNode#setParent(AtomTreeNodeGroup)
	 */
	public void setParent(AtomTreeNodeGroup parent) {
		super.setParent(parent);
		childIterator.reset();
		while(childIterator.hasNext()) childIterator.nextAtom().node.setParent(this);//notifies children that there was a change in the hierarchy
	}            
	
    public final Atom firstChildAtom() {return childList.getFirst();}
    public final Atom lastChildAtom() {return childList.getLast();}

    public Atom randomAtom() {return childList.getRandom();}
    
    /**
     * Indicates whether the children of this group are themselves atom groups,
     * or are leaf atoms.  Returns true if any child atom is a group, and thus does
     * not imply that all child atoms are groups.
     */
    public boolean childrenAreGroups() {
        return ((AtomTypeGroup)atom.type).childrenAreGroups;
   //     return !firstChildAtom().node.isLeaf();
    }

    /**
     * Gets the child atom corresponding to the given index, numbering the first atom as zero
     * and the last atom as Count-1.
     */
    public Atom getAtom(int index) {
        return childList.get(index);
    }//end of getAtom
            

    public boolean isLeaf() {return false;}
    
    public Atom firstLeafAtom() {
        return findFirstLeafAtom();//firstLeafAtom;  //firstLeafAtom is not updated correctly
                                                     //fix this if using stored value
    }
    private Atom findFirstLeafAtom() {
        if(childrenAreGroups()) {
            childIterator.reset();
            while(childIterator.hasNext()) {
                Atom a1 = childIterator.nextAtom().node.firstLeafAtom();
                if(a1 != null) return a1;
            }
            return null;
        }
          /* using iterator for loop seems to cause problem with return of null in AtomPairIterator inner loop even though it reports hasNext=true  
            childIterator.reset(UP);
            while(childIterator.hasNext()) {
                Atom atom = ((AtomGroup)childIterator.next()).firstLeafAtom();
                if(atom != null) return atom;
            }
            return null;
        }*/
        else return firstChildAtom();
    }
    
    /**
     * Returns the last leaf atom descended from this group.
     */
    public Atom lastLeafAtom() {
        if(childrenAreGroups()) {
            AtomLinker a0 = childList.header.previous;
            while(a0 != childList.header) {
                if(a0.atom == null) {//skip tabs
                    a0 = a0.previous;
                } else {
                    Atom a1 = a0.atom.node.lastLeafAtom();
                    if(a1 != null) return a1;
                    else a0 = a0.previous;
                }
            }
            return null;
        }
        else return lastChildAtom();
    }
    
    /**
     * Rearranges the hierarchy of atoms below this node, resulting
     * in a situation where all children of this node are leaf atoms.
     */
     //need to find way to update childrenAreGroups to give false 
     //after this method is invoked
     
     //commented because this is just a detailed sketch of the method; not tested at all
    /*
    public void flatten() {
        AtomIteratorListSimple iterator = new AtomIteratorListSimple();
        AtomList tempChildList = new AtomList(childList);//copy because we're modifying the child list
        childList.clear();
        childIterator.setBasis(tempChildList);
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom child = childIterator.next();
            if(child.node instanceof AtomTreeGroupNode) {
                AtomTreeNodeGroup childNode = (AtomTreeNodeGroup)child.node;
                childNode.flatten();
                iterator.setBasis(childNode.childList);
                iterator.reset();
                while(iterator.hasNext()) {
                    childList.add(iterator.next().seq);
                }
            } else {
                childList.add(child.seq)
            }
        }
    }//end of flatten
    
    */    
    public int leafAtomCount() {return leafAtomCount;}
    public int childAtomCount() {return childList.size();}

    /**
     * Returns the children of this group in an array of atoms.
     * Array is constructed on-the-fly, and is not updated with any
     * subsequent atom addition/removals.  Since array construction is
     * involved, this method can be expensive in computationally intensive
     * situations involving repeated calls (this should be avoided).
     */
    public Atom[] childAtomArray() {
        return childList.toArray();
    }
    
    public void removeAllChildren() {
        Atom[] array = childAtomArray();
        for(int i=0; i<array.length; i++) {
            array[i].sendToReservoir();
        }
    }
    
    public boolean isResizable() {return true;}
    
    /**
     * Notifies this atom group that an atom has been added to it or one of its descendants.
     */
    public void addAtomNotify(Atom atom) {
        leafAtomCount += atom.node.leafAtomCount();
        parentNode().addAtomNotify(atom);
    }
    
    public void removeAtomNotify(Atom atom) {
        leafAtomCount -= atom.node.leafAtomCount();
        parentNode().removeAtomNotify(atom);
    }
      
    protected int leafAtomCount;
    
    //childlist is public, but should not add/remove atoms except via node's methods.
    //consider a mechanism to ensure this; a inner mutator class made available only
    //to list's creator, for example (still wouldn't prevent modification via direct
    //access of entry classes).
    public final AtomList childList = new AtomList();
    private final AtomIteratorListSimple childIterator = new AtomIteratorListSimple(childList);
        
    public static final AtomTreeNode.Factory FACTORY = new AtomTreeNode.Factory() {
        public AtomTreeNode makeNode(Atom atom, AtomTreeNodeGroup parent) {
            return new AtomTreeNodeGroup(atom, parent);
        }
    };
    
}