/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.AtomTypeGroup;
import etomica.AtomType.SphericalTop;


/**
 * Atom type for a sphere that has some feature depending upon an orientation coordinate.
 * For example an orientational dependent potential may be attached to an otherwise spherical atom
 */
public final class AtomTypeOrientedSphere extends AtomTypeSphere implements SphericalTop {
    
    private final double[] I = new double[3];
    public AtomTypeOrientedSphere(AtomTypeGroup parentType, double m, double d) {
        super(parentType, m, d);
        updateI();
    }
    public double[] momentOfInertia() {return I;}
    
    private void updateI() {
        if(I == null) return;
        I[0] = 0.4*this.getMass()*radius*radius;  //moment of inertia of a sphere = 2/5 m R^2 (should modify to arbitrary dimension)
        I[1] = I[0];
        I[2] = I[1];
    }
    
    public void setMass(double m) {
        super.setMass(m);
        updateI();
    }
    public void setDiameter(double d) {
        super.setDiameter(d);
        updateI();
    }
}