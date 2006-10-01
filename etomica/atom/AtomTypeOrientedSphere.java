/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.atom.AtomType.SphericalTop;
import etomica.chem.elements.Element;


/**
 * Atom type for a sphere that has some feature depending upon an orientation coordinate.
 * For example an orientational dependent potential may be attached to an otherwise spherical atom
 */
public final class AtomTypeOrientedSphere extends AtomTypeSphere implements SphericalTop {
    
    private final double[] I = new double[3];
    public AtomTypeOrientedSphere(Element element, double d) {
        super(element, d);
        updateI();
    }
    public double[] momentOfInertia() {return I;}
    
    private void updateI() {
        if(I == null) return;
        I[0] = 0.1*this.getMass()*diameter*diameter;  //moment of inertia of a sphere = 2/5 m R^2 (should modify to arbitrary dimension)
        I[1] = I[0];
        I[2] = I[1];
    }
    
    public void setDiameter(double d) {
        super.setDiameter(d);
        updateI();
    }
}