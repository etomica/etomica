package etomica.action;

import etomica.api.IAtomList;

/**
 * Interface for a class that can perform an action on an atom set.
 */
public interface AtomsetAction {
    
	public void actionPerformed(IAtomList atoms);
}
