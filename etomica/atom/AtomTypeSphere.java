/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.util.Default;

/**
 * Atom type for a simple monatomic atom that has a length scale associated
 * with its size.  Position definition is the atom's coordinate 
 * (AtomPositionDefinitionSimple).
 */

public class AtomTypeSphere extends AtomTypeLeaf {
    
    double diameter, radius;
    
    public AtomTypeSphere(AtomTypeGroup parentType) {
        this(parentType, Default.ATOM_MASS, Default.ATOM_SIZE);
    }
    public AtomTypeSphere(AtomTypeGroup parentType, double m, double d) {
        super(parentType, m);
        setDiameter(d);
    }
                
    public double diameter(Atom a) {return diameter;}
    public double radius(Atom a) {return radius;}
    
    /**
    * Sets diameter of this atom and updates radius accordingly.
    *
    * @param d   new value for diameter
    */
    public void setDiameter(double d) {diameter = d; radius = 0.5*d;}
}