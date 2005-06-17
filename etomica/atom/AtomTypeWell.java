/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.AtomTypeGroup;

public final class AtomTypeWell extends AtomTypeSphere {  
    
    private double lambda;                    //diameter of well, in units of core diameter
    private double wellDiameter, wellRadius;  //size of well, in simulation units
    
    public AtomTypeWell(AtomTypeGroup parentType, double m, double d, double l) {
        super(parentType, m, d);
        setDiameter(d);
        setLambda(l);
    }
                
    public final double lambda() {return lambda;}
    public final double wellDiameter() {return wellDiameter;}
    public final double wellRadius() {return wellRadius;}
    
    public final void setDiameter(double d) {
        super.setDiameter(d); 
        setLambda(lambda);
    }
    public final void setLambda(double l) {
        lambda = l; 
        wellDiameter = lambda*diameter; 
        wellRadius = 0.5*wellDiameter;
    }
}