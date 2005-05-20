/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.Atom;
import etomica.AtomType;
import etomica.Default;
import etomica.Parameter;


public class AtomTypeSphereVariable extends AtomTypeSphere {
    
    int index;
    
    public AtomTypeSphereVariable(AtomType parentType, int index) {
        this(parentType, index, Default.ATOM_MASS);
    }
    public AtomTypeSphereVariable(AtomType parentType, int index, double m) {
        super(parentType, m, Double.NaN);
        this.index = index;
    }
                
    public double diameter(Atom a) {return ((Parameter.Size)(a.allatomAgents[index])).getSigma();}
    public double radius(Atom a) {return 0.5*diameter(a);}
    
    /**
    * Sets diameter of this atom and updates radius accordingly.
    *
    * @param d   new value for diameter
    */
    public void setDiameter(double d) {if(!Double.isNaN(d)) throw new RuntimeException("Unexpected call to AtomType.SphereVariable.setDiameter method");}
}