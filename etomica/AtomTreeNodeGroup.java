package etomica;

/**
 * Implementation of AtomTreeNode for non-leaf node.
 */

public class AtomTreeNodeGroup implements AtomTreeNode {
    
    public AtomTreeNodeGroup(Atom atom) {
        this.atom = atom;
    }

    public Atom atom() {return atom;}
        
    public AtomGroup parentGroup() {
        return parentGroup;
    }
    public AtomTreeNodeGroup parentNode() {
        return parentNode;
    }
    
    public Class childSequencerClass() {
        try { return firstChildAtom().seq.getClass();}
        catch(NullPointerException e) {return null;}//in case firstChild is null
    }
    
    /**
     * Returns the molecule in which this atom resides.  A "molecule" is an atomgroup
     * that is one step below a species agent in the hierarchy of atomgroups.
     */
    public Atom parentMolecule() {
        return (parentNode.atom() instanceof SpeciesAgent) ? this.atom : parentNode.parentMolecule();
    }
    
    public void setParentGroup(AtomGroup parent) {
        parentGroup = parent;
        parentNode = (parent != null) ? (AtomTreeNodeGroup)parent.node : null;
        if(parentNode != null) {
            depth = parentNode.depth() + 1;
            parentPhase = parentNode.parentPhase();
        }
    }

    public void setDepth(int d) {
        depth = d;
        for(Atom atom=firstChildAtom(); atom!=null; atom=atom.nextAtom()) {
            atom.node.setDepth(d+1);
            if(atom == lastChildAtom()) break;
        }
        if(parentNode != null) parentPhase = parentNode.parentPhase();
    }
    
    public final Atom firstChildAtom() {return childList.getFirst();}
    public final Atom lastChildAtom() {return childList.getLast();}

    public Atom randomAtom() {return childList.getRandom();}
    
    /**
     * Indicates whether the children of this group are themselves atom groups,
     * or are leaf atoms.
     */
    public boolean childrenAreGroups() {
        return !firstChildAtom().node.isLeaf();
    }

    /**
     * Gets the child atom corresponding to the given index, numbering the first atom as zero
     * and the last atom as Count-1.
     */
    public Atom getAtom(int index) {
        return childList.get(index);
    }//end of getAtom
            
    /**
     * Simulation in which this atom resides
     */
    public Simulation parentSimulation() {return parentPhase().parentSimulation();}        
    /**
     * Phase in which this atom resides
     */
    public Phase parentPhase() {return parentPhase;}//parentNode.parentPhase();}

    public Species parentSpecies() {return parentNode.parentSpecies();}
    
    public SpeciesAgent parentSpeciesAgent() {return parentNode.parentSpeciesAgent();}

    /**
     * Returns the depth of this atom in the atom hierarchy.  That is, returns
     * the number of parent relations between this atom and the species master.
     */
    public int depth() {return depth;}//return (parentGroup != null) ? parentGroup.depth()+1 : 0;}
    
    public boolean isLeaf() {return false;}
    
    /**
     * Integer assigned to this atom by its parent molecule.
     * Assigned during construction of atom.
     */
    public final int index() {return atomIndex;}
    public final void setIndex(int i) {atomIndex = i;}

    /**
     * Returns true if the given atom is in the hierarchy of parents of this atom,
     * or if the given atom is this atom.  Returns true, for example, if the given
     * atom is this atom's parent, or its parent's parent, etc.
     */ 
    public boolean isDescendedFrom(Atom group) {
        return (this.atom == group) || (parentNode != null && parentNode.isDescendedFrom(group));
    }
     /*   AtomGroup ancestor = parentGroup;
        while(ancestor != null) {
            if(ancestor == group) return true;
            ancestor = ancestor.parentGroup();
        }
        return false;
    }*/
        
    public Atom firstLeafAtom() {
        return findFirstLeafAtom();//firstLeafAtom;  //firstLeafAtom is not updated correctly
                                                     //fix this if using stored value
    }
    private Atom findFirstLeafAtom() {
        if(childrenAreGroups()) {
            for(Atom a0=firstChildAtom(); a0!=null; a0=a0.nextAtom()) {
                Atom a1 = a0.node.firstLeafAtom();
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
            for(Atom a0=lastChildAtom(); a0!=null; a0=a0.previousAtom()) {
                Atom a1 = ((AtomGroup)a0).node.lastLeafAtom();
                if(a1 != null) return a1;
            }
            return null;
        }
          /* using iterator for loop seems to cause problem with return of null in AtomPairIterator inner loop even though it reports hasNext=true  
            childIterator.reset(DOWN);
            while(childIterator.hasNext()) {
                Atom atom = ((AtomGroup)childIterator.next()).lastLeafAtom();
                if(atom != null) return atom;
            }
            return null;
        } */
        else return lastChildAtom();
    }
    
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
    
            
    public void addAtom(Atom aNew) {
        if(!resizable) return; //should define an exception
        //ensure that the given atom is compatible
        //creator() is not defined for speciesAgent, so this creates a problem when adding
        //new molecule; commented out vetoAddition call until this is fixed or abandoned.
        if(aNew == null /* || creator().vetoAddition(aNew)*/) { //should ensure group/atom consistency with other children
            throw new NullPointerException("Attempting to add null Atom in AtomTreeNodeGroup.addAtom");
        }
        //set up links within this group
        aNew.node.setParentGroup((AtomGroup)atom);
        aNew.node.setDepth(depth+1);
        
        if(childAtomCount() > 0) { //siblings
            aNew.node.setIndex(lastChildAtom().node.index()+1);
        }
        else {  //only child
            aNew.node.setIndex(0);
        }  
        
        childList.add(aNew.seq);
        
        //set up leaf links
/*        if(aNew instanceof AtomGroup) {
            AtomGroup aNewGroup = (AtomGroup)aNew;//copy to save multiple casts in calls to AtomGroup methods
            Atom lastNewLeaf = aNewGroup.node.lastLeafAtom();
            if(lastNewLeaf != null) {//also implies firstNewLeaf != null
                lastNewLeaf.seq.setNextAtom(findNextLeafAtom(aNewGroup));
                Atom previous = findPreviousLeafAtom(aNewGroup);
                if(previous != null) previous.seq.setNextAtom(aNewGroup.node.firstLeafAtom());
            } //no leaf changes anywhere if aNew has no leaf atoms 
        }
        else {//aNew is a leaf, not a group
                //set next in case it wasn't set by aNew.setNextAtom(lastChildAtom.next()) line, above
            if(aNew.nextAtom() == null) aNew.seq.setNextAtom(findNextLeafAtom(aNew.node.parentGroup()));
               //set previous in case it wasn't set by lastChildAtom.setNextAtom(aNew) line, above
            if(aNew.previousAtom() == null) {
                Atom previous = findPreviousLeafAtom(aNew.node.parentGroup());
                if(previous != null) previous.seq.setNextAtom(aNew);
            }
        }
        */
        addAtomNotify(aNew);
    }//end of addAtom


    public void removeAtom(Atom a) {
        if(a == null || !resizable) return;
        //ensure that the given atom is compatible
        if(a.node.parentGroup() != atom) {
            System.out.println("Error in AtomGroup.removeAtom:  See source code");  //should throw an exception
            return;
        }
        childList.remove(a);
        //childList.remove(a.seq); //could do this if method not private in AtomList
        
        removeAtomNotify(a);
        a.node.setParentGroup(null);//must follow notify for SpeciesMaster moleculeCount to be updated
    }//end of removeAtom

    /**
     * Removes all child atoms of this group.  Maintains their
     * internal links, and returns the first child, so they may be recovered 
     * for use elsewhere.
     * Does not remove this group from its parent group.
     */
     //doesn't perform removeAtomNotify
    public Atom removeAll() {
        //save first child for return, then erase links to them
        Atom first = firstChildAtom();
        childList.clear();
        return first;
    }//end of removeAll

    /**
     * Notifies this atom group that an atom has been added to it or one of its descendants.
     */
    public void addAtomNotify(Atom atom) {
        leafAtomCount += atom.node.leafAtomCount();
        if(parentGroup() != null) parentGroup().node.addAtomNotify(atom);
    }
    
    public void removeAtomNotify(Atom atom) {
        leafAtomCount -= atom.node.leafAtomCount();
        if(parentGroup() != null) parentGroup().node.removeAtomNotify(atom);
    }
    

    /**
     * Searches the atom hierarchy up from (and not including) the given group 
     * to find the closest (leaf) atom in that direction.
     */
    private static Atom findNextLeafAtom(AtomGroup a) {
        AtomGroup parent = a; //search up siblings first
        while(parent != null) {
            AtomGroup uncle = (AtomGroup)parent.nextAtom();
            while(uncle != null) {
                Atom first = uncle.node.firstLeafAtom();
                if(first != null) return first; //found it
                uncle = (AtomGroup)uncle.nextAtom();
            }
            parent = parent.node.parentGroup();
        }
        return null;
    }//end of findNextLeafAtom
    
    /**
     * Searches the atom hierarchy down from (and not including) the given group 
     * to find the closest (leaf) atom in that direction.
     * The first leaf atom of the given group will have the returned atom as its previous atom. 
     */
    private static Atom findPreviousLeafAtom(AtomGroup a) {
        AtomGroup parent = a;
        while(parent != null) {
            AtomGroup uncle = (AtomGroup)parent.previousAtom();
            while(uncle != null) {
                Atom last = uncle.node.lastLeafAtom();
                if(last != null) return last;
                uncle = (AtomGroup)uncle.previousAtom();
            }
            parent = parent.node.parentGroup();
        }
        return null;
    }//end of findPreviousLeafAtom

    protected final Atom atom;
    protected int depth;
    protected int atomIndex;
    protected int leafAtomCount;
    private AtomTreeNodeGroup parentNode;
    private AtomGroup parentGroup;
    private Phase parentPhase;
    
    //childlist is public, but should not add/remove atoms except via node's methods.
    //consider a mechanism to ensure this; a inner mutator class made available only
    //to list's creator, for example (still wouldn't prevent modification via direct
    //access of entry classes).
    public final AtomList childList = new AtomList();
    private boolean resizable = true;
        
    public static final AtomTreeNode.Factory FACTORY = new AtomTreeNode.Factory() {
        public AtomTreeNode makeNode(Atom atom) {
            return new AtomTreeNodeGroup(atom);
        }
    };
    
}