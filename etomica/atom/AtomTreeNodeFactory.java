package etomica.atom;


/**
 * Interface for class that makes an AtomTreeNode.  Classes implementing this
 * interface are given to the Atom constructor to make the tree-node for the
 * atom.  Most implementations of this interface are simple inner classes defined
 * in the AtomTreeNode subclass that they construct.
 *
 * @author David Kofke
 *
 */
/*
 * History
 * Created on Feb 12, 2005 by kofke
 */
public interface AtomTreeNodeFactory {
    AtomTreeNode makeNode(Atom atom);
}