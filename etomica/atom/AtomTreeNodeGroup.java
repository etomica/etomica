package etomica.atom;

import etomica.Atom;

/**
 * Implementation of AtomTreeNode for non-leaf node.
 */
  
public class AtomTreeNodeGroup extends AtomTreeNode {
    
    public AtomTreeNodeGroup(Atom atom) {
        super(atom);
    }
        
    public int newChildIndex() {
        if(!childList.isEmpty()) { //siblings
            return (childList.getLast().node.index()+1);
        }
        else {  //only child
            return(0);
        }
    }

	/**
	 * Invoke superclass method and then notify children.
	 * 
	 * @see etomica.atom.AtomTreeNode#setParent(AtomTreeNodeGroup)
	 */
	public void setParent(AtomTreeNodeGroup parent) {
		super.setParent(parent);
        AtomLinker link = childList.header.next;
        while(link != childList.header) {
            link.atom.node.setParent(this);//notifies children that there was a change in the hierarchy
            link = link.next;
        }
    }
	
    public boolean isLeaf() {return false;}
    
    public Atom firstLeafAtom() {
        AtomLinker link = childList.header.next;
        while(link != childList.header) {
            Atom a1 = link.atom.node.firstLeafAtom();
            if(a1 != null) return a1;
            else link = link.next;
        }
        return null;
    }
    
    /**
     * Returns the last leaf atom descended from this group.
     */
    public Atom lastLeafAtom() {
        AtomLinker link = childList.header.previous;
        while(link != childList.header) {
            Atom a1 = link.atom.node.lastLeafAtom();
            if(a1 != null) return a1;
            else link = link.previous;
        }
        return null;
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
            array[i].node.dispose();
        }
    }
    
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
        
    public static final AtomTreeNodeFactory FACTORY = new AtomTreeNodeFactory() {
        public AtomTreeNode makeNode(Atom atom) {
            return new AtomTreeNodeGroup(atom);
        }
    };
    
}