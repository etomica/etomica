package etomica.normalmode;

import etomica.api.IBox;

/**
 * Preliminary interface for proposed IPotential interface
 */
public interface IAPIPotential {

    /**
     * Returns the potential energy of the given box.
     */
    public double calculateEnergy(IBox box);
}
