package etomica.action;

import etomica.atom.AtomSet;
import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;

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
 
public abstract class AtomActionAdapter extends AtomsetActionAdapter implements AtomAction {
    
    public AtomActionAdapter() {
        atoms = new AtomSetSinglet();
    }
    
    public void setAtom(IAtom a) {((AtomSetSinglet)atoms).atom = a;}
    public IAtom getAtom() {return ((AtomSetSinglet)atoms).atom;}

    public void actionPerformed(AtomSet a) {
        actionPerformed(a.getAtom(0));
    }
    
    /**
     * Method that defines the action to be performed on the atom
     * @param a Atom passed to method by iterator
     */
    public abstract void actionPerformed(IAtom a);
    
}   