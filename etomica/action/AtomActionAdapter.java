package etomica.action;

import etomica.Atom;

/**
 * Base class for classes that perform some elementary action on an atom.
 * These classes see use in several ways:
 * <ul>
 * <li>they can be passed to the allAtoms method of atom iterators, which then performs
 * the specified action on all the atoms of the iterator.  
 * <li>they can be used in a DisplayPhase.AtomActionWrapper, to specify some action in response to 
 * selection of an atom by the mouse.  
 * <li>they may be used to generate a Monte Carlo trial in an MCMove object.
 * 
 * @author David Kofke
 * @see Atom.Iterator
 * @see DisplayPhase.AtomActionWrapper
 */
 
public abstract class AtomActionAdapter extends AtomsetActionAdapter implements AtomAction {
    
    public void setAtom(Atom a) {atoms[0] = a;}
    public Atom getAtom() {return atoms[0];}

    public void actionPerformed(Atom[] a) {
        actionPerformed(a[0]);
    }
    
    /**
     * Method that defines the action to be performed on the atom
     * @param a Atom passed to method by iterator
     */
    public abstract void actionPerformed(Atom a);
    
}   