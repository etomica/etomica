package etomica.atom;

import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.chem.elements.Element;
import etomica.space.ISpace;


/**
 * Atom type for a sphere that has some feature depending upon an orientation coordinate.
 * For example an orientational dependent potential may be attached to an otherwise spherical atom
 */
public class AtomTypeOrientedSphere extends AtomTypeSphere implements IAtomTypeOriented {
    
    protected final IVectorMutable I;
    public AtomTypeOrientedSphere(Element element, double d, ISpace space) {
        super(element, d);
        I = space.makeVector();
        updateI();
    }
    public IVector getMomentOfInertia() {return I;}
    
    protected void updateI() {
        //moment of inertia of a sphere = 2/5 m R^2 (should modify to arbitrary dimension)
        if(I != null)
            I.E(0.1*this.getMass()*diameter*diameter);
    }
    
    public void setDiameter(double d) {
        super.setDiameter(d);
        updateI();
    }
}