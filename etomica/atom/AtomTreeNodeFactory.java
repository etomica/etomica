package etomica.atom;

import etomica.Atom;



/*
 * History
 * Created on Feb 12, 2005 by kofke
 */
public interface AtomTreeNodeFactory {
    public AtomTreeNode makeNode(Atom atom);
}