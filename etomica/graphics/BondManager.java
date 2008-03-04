package etomica.graphics;

import etomica.api.IAtomSet;

/**
 * The BondManager is responsible for creation and disposal of bonds between
 * atom pairs.
 *
 * @author Andrew Schultz
 */
public interface BondManager {
    
    /**
     * Creates a bond object for the given pair with type determined by the
     * bondType object.  The returned object can be used later to
     * tell the BondManager to remove the bond.
     */
    public Object makeBond(IAtomSet pair, Object bondType);
    
    /**
     * Notifies the BondManager that the given bond no longer exists in the
     * system.  The given bond object must have been returned from a previous
     * call to makeBond.
     */
    public void releaseBond(Object bond);
}