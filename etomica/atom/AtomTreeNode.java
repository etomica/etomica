package etomica.atom;

import etomica.Atom;
import etomica.Phase;
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
 
public abstract class AtomTreeNode {
    
    public AtomTreeNode(Atom atom) {
        this.atom = atom;
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
    	if(parent != null && parent == parentNode) {
			depth = parentNode.depth() + 1;
			return;
        }
    	
        //old parent is not null; remove from it
        if(parentNode != null) {
            parentNode.childList.remove(atom.seq);        
            parentNode.removeAtomNotify(atom);
        }
        
        parentNode = parent;

        if(parentNode == null) {//new parent is null
            return;
        }

        //new parent is not null
        depth = parentNode.depth() + 1;

        setIndex(parentNode.newChildIndex());
        
        parentNode.childList.add(atom.seq);
        parentNode.addAtomNotify(atom);
    }//end of addAtom

    public Atom atom() {return atom;}
        
    public Atom parentGroup() {
        return (parentNode != null) ? parentNode.atom : null;
    }
    
    public AtomTreeNodeGroup parentNode() {
        return parentNode;
    }

    /**
     * Returns the molecule in which this atom resides.  A "molecule" is an atomgroup
     * that is one step below a species agent in the hierarchy of atomgroups.
     */
    public Atom parentMolecule() {
        if(parentNode == null) return null;
        return (parentNode instanceof SpeciesAgent.AgentAtomTreeNode) ? this.atom : parentNode.parentMolecule();
    }
                
    /**
     * Phase in which this atom resides
     */
    public Phase parentPhase() {
        return (parentNode != null) ? parentNode.parentPhase() : null;
    }
    
    public SpeciesAgent parentSpeciesAgent() {
        return (parentNode != null) ? parentNode.parentSpeciesAgent() : null;
    }

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
        if(parentNode == null) return null;
        return (parentNode == node) ? this : parentNode.childWhereDescendedFrom(node);
    }
    
    
    /**
     * Returns true if this atom preceeds the given atom in the atom sequence.
     * Returns false if the given atom is this atom, or (of course) if the
     * given atom instead preceeds this one.
     */
    public boolean preceeds(AtomTreeNode node) {
        //want to return false if atoms are the same atoms
        if(parentNode == null) return false;
        if(node == null) return true;
        if(this.parentNode == node.parentNode) return atomIndex < node.atomIndex;//works also if both parentGroups are null
        if(depth == node.depth) return parentNode.preceeds(node.parentNode);
        if(depth < node.depth) return this.preceeds(node.parentNode);
        /*if(this.depth > atom.depth)*/ return parentNode.preceeds(node);
    }
    
    public boolean preceeds(Atom a) {
        return (a != null) ? preceeds(a.node) : true;
    }
        
    protected final Atom atom;
    protected int depth;
    protected int atomIndex;
    private AtomTreeNodeGroup parentNode;
    
}//end of AtomTreeNode