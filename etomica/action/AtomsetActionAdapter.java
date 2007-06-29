package etomica.action;

import etomica.atom.AtomSet;

/**
 * Base class for classes that perform some elementary action on an atom.
 * These classes see use in several ways:
 * <ul>
 * <li>they can be passed to the allAtoms method of atom iterators, which then performs
 * the specified action on all the atoms of the iterator.  
 * <li>they can be used in a DisplayBox.AtomActionWrapper, to specify some action in response to 
 * selection of an atom by the mouse.  
 * <li>they may be used to generate a Monte Carlo trial in an MCMove object.
 * 
 * @author David Kofke
 */
 
public abstract class AtomsetActionAdapter implements AtomsetAction, java.io.Serializable {
    
    protected AtomSet atoms;

    public void setAtoms(AtomSet a) {atoms = a;}
    public AtomSet getAtoms() {return atoms;}

    /**
     * Performs the defined action using the atom most recently specified by setAtom or by the last call to actionPerformed(Atom a).
     * Performs no action if the atom is null.
     */
    public void actionPerformed() {
        if(atoms != null) actionPerformed(atoms);
    }
    
    /**
     * Method that defines the action to be performed on the atom
     * @param a Atom passed to method by iterator
     */
    public abstract void actionPerformed(AtomSet a);
}   
