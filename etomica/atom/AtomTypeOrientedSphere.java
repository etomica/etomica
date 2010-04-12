package etomica.atom;

import etomica.api.IElement;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;


/**
 * Atom type for a sphere that has some feature depending upon an orientation coordinate.
 * For example an orientational dependent potential may be attached to an otherwise spherical atom
 */
public class AtomTypeOrientedSphere extends AtomTypeLeaf implements IAtomTypeOriented {
    
    protected final IVectorMutable I;
    public AtomTypeOrientedSphere(IElement element, ISpace space) {
        super(element);
        I = space.makeVector();
    }
    public IVector getMomentOfInertia() {return I;}

    public void setMomentOfInertia(double moment) {
        I.E(moment);
    }
}