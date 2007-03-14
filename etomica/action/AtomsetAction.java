package etomica.action;

import etomica.atom.AtomSet;

/**
 * Interface for a class that can perform an action on an atom set.
 */
public interface AtomsetAction extends Action {
    
	public void actionPerformed(AtomSet atoms);
    public void setAtoms(AtomSet atom);
    public AtomSet getAtoms();
    
}
