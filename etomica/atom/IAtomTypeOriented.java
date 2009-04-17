package etomica.atom;

import etomica.api.IAtomType;
import etomica.api.IVector;

public interface IAtomTypeOriented extends IAtomType {

    /**
     * Returns the principle components of the moment of inertia of the
     * atom within the body-fixed frame.  Do NOT modify the returned moment
     * of inertia returned.
     */
    public IVector getMomentOfInertia();
}
