/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica;

import etomica.Parameter.Size;


public class AtomTypeSphereVariable extends AtomTypeSphere {
    
    int index;
    
    public AtomTypeSphereVariable(AtomFactory creator, int index) {
        this(creator, index, Default.ATOM_MASS);
    }
    public AtomTypeSphereVariable(AtomFactory creator, int index, double m) {
        super(creator, m, Double.NaN);
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