package etomica.atom;

import etomica.api.IMolecule;
import etomica.api.IVectorMutable;


/**
 * Interface for an molecule that holds vectors for velocity.
 */
public interface IMoleculeKinetic extends IMolecule {

    public IVectorMutable getVelocity();
 
}
