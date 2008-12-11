package etomica.junit.atom;

import etomica.api.IMoleculeList;

/**
 * Interface for a class that can perform an action on an atom set.
 */
public interface MoleculesetAction {
    
	public void actionPerformed(IMoleculeList atoms);
}
