package etomica.action;

import etomica.api.IMolecule;


/**
 * Interface for a class that can perform an action on an atom.
 */
public interface MoleculeAction {

    public void actionPerformed(IMolecule atom);
}
