package etomica.atom;

import etomica.Atom;
import etomica.Phase;
import etomica.Species;
import etomica.SpeciesAgent;

/**
 * Interface for a node in the atom tree.  Atoms are arranged in a tree structure 
 * to define collections of molecules, individual molecules, atom groups and atoms.  
 * All of these elements are represented by the Atom class, and the structure that
 * forms the larger elements from the smaller ones is maintained by the node class
 * associated with each Atom.  The node for each atom holds information about its
 * parent and its children (if any).  Leaf nodes have no children, but they are defined
 * to identify themselves as their only child (this design makes certain types of
 * pair loops easier to perform).<br>
 *
 * (The node for) a SpeciesMaster instance is the root of the atom tree.  Directly below
 * it are instances of SpeciesAgent class, which are constructed each Phase by each Species.
 * Below the SpeciesAgent are instances of Atoms that are logically the molecules of the
 * system.  The molecules may be leaf atoms, ending the tree, or they may contain child 
 * atoms (or other groups of atoms) used to form the molecule.
 *
 * @author David Kofke
 */
 
 /* History
  * 11/15/02 (DAK/DW) modified parentPhase method
  * 12/04/02 (DAK) added isDescendedFrom(AtomTreeNodeGroup) method, and modified
  *           isDescendedFrom(Atom) method to work through it.
  * 12/06/02 (DAK) added childWhereDescendedFrom method.
  * 08/12/03 (DAK) modifed as indicated in comments in code; modified dispose()
  * method to operate via setParent method
  */
 

public abstract class AtomTreeNode {
    
    public AtomTreeNode(Atom atom, AtomTreeNodeGroup parent) {
        this.atom = atom;
//        if(parent == null && !(atom instanceof SpeciesMaster || atom instanceof AtomReservoir)) {//only SpeciesMaster and AtomReservoir can have null parent
//            throw new NullPointerException("Illegal specification of null parent in AtomTreeNode constructor");
//        } else if(parent != null) {
        	if(parent != null) {
            //(add check that parent is resizable)
            parentNode = parent;
            parentGroup = parent.atom;
            depth = parent.depth() + 1;
            parentPhase = parent.parentPhase();
            setIndex(parent.newChildIndex());
            parent.childList.add(atom.seq);
 //           parent.addAtomNotify(atom); //invoked instead in constructor of atom in which this node is being placed
        }
    }
    
    public abstract boolean isLeaf();
    public abstract Atom firstLeafAtom();
    
    /**
     * Returns the last leaf atom descended from this group.
     */
    public abstract Atom lastLeafAtom();
    
    public abstract int leafAtomCount();
    public abstract int childAtomCount();
    
    public void dispose() {
    	setParent((AtomTreeNodeGroup)null);
    }

    public void setParent(Atom parent) {
        setParent(parent==null ? (AtomTreeNodeGroup)null : (AtomTreeNodeGroup)parent.node);
    }
    
    public void setParent(AtomTreeNodeGroup parent) {
    	
    	//parent isn't changing, but may need to update fields (added this 'if' block 08/12/03 (DAK))
    	if(parent == parentNode) {
			parentPhase = parentNode.parentPhase();
			depth = parentNode.depth() + 1;
			return;
        }
    	
        if(parent != null && !parent.isResizable()) return;  //throw an exception?
        
        //old parent is not null; remove from it
        if(parentNode != null) {
            if(!parentNode.isResizable()) return;//exception
            parentNode.childList.remove(atom.seq);        
            parentNode.removeAtomNotify(atom);
        }
        
        parentNode = parent;

        if(parentNode == null) {//new parent is null
            parentGroup = null;
            parentPhase = null;
            return;
        }

        //new parent is not null
        parentGroup = parentNode.atom;
        parentPhase = parentNode.parentPhase();
        depth = parentNode.depth() + 1;

        setIndex(parentNode.newChildIndex());
        
        parentNode.childList.add(atom.seq);
        
        //should notify this node's children of change

        parentNode.addAtomNotify(atom);
    }//end of addAtom

    public Atom atom() {return atom;}
        
    public Atom parentGroup() {
        return parentGroup;
    }
    public AtomTreeNodeGroup parentNode() {
        return parentNode;
    }

    /**
     * Returns the molecule in which this atom resides.  A "molecule" is an atomgroup
     * that is one step below a species agent in the hierarchy of atomgroups.
     */
    public Atom parentMolecule() {
        return (parentNode.atom() instanceof SpeciesAgent) ? this.atom : parentNode.parentMolecule();
    }
                
    /**
     * Phase in which this atom resides
     */
    public Phase parentPhase() {//return parentPhase;}//parentNode.parentPhase();}
        return (parentPhase != null) ? parentPhase : 
                    (parentNode != null) ? parentNode.parentPhase() : null;
    }

    public Species parentSpecies() {return parentNode.parentSpecies();}
    
    public SpeciesAgent parentSpeciesAgent() {return parentNode.parentSpeciesAgent();}

    /**
     * Returns the depth of this atom in the atom hierarchy.  That is, returns
     * the number of parent relations between this atom and the species master.
     */
    public final int depth() {return depth;}//return (parentGroup != null) ? parentGroup.depth()+1 : 0;}
    
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
    	if(group == null) return false;
        else return this.isDescendedFrom(group.node);
     //   return (this.atom == group) || (parentNode != null && parentNode.isDescendedFrom(group));
    }
    
    public boolean isDescendedFrom(AtomTreeNode node) {
    	if(node == null) return false;
        else return (this == node) || (parentNode != null && parentNode.isDescendedFrom(node));
    }
    
    /**
     * Returns the child of the given node through which this node is derived.
     * If given node is parent node of this, returns this.
     * If this node is not descended from the given node, returns null.
     */
    public AtomTreeNode childWhereDescendedFrom(AtomTreeNode node) {
        return (parentNode == node) ? this : parentNode.childWhereDescendedFrom(node);
    }
    
    
    /**
     * Returns true if this atom preceeds the given atom in the atom sequence.
     * Returns false if the given atom is this atom, or (of course) if the
     * given atom instead preceeds this one.
     */
    public boolean preceeds(Atom a) {
        //want to return false if atoms are the same atoms
        if(a == null) return true;
        if(parentGroup == a.node.parentGroup) return atomIndex < a.node.atomIndex;//works also if both parentGroups are null
        if(depth == a.node.depth) return parentNode.preceeds(a.node.parentGroup);
        if(depth < a.node.depth) return this.preceeds(a.node.parentGroup);
        /*if(this.depth > atom.depth)*/ return parentNode.preceeds(a);
    }
        
    protected final Atom atom;
    protected int depth;
    protected int atomIndex;
    private AtomTreeNodeGroup parentNode;
    private Atom parentGroup;
    private Phase parentPhase;

    
    public interface Factory {
        public AtomTreeNode makeNode(Atom atom, AtomTreeNodeGroup parent);
    }
    
}//end of AtomTreeNode