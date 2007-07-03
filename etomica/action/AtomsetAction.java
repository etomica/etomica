package etomica.action;

import etomica.atom.AtomSet;

/**
 * Interface for a class that can perform an action on an atom set.
 */
public interface AtomsetAction {
    
	public void actionPerformed(AtomSet atoms);
}
