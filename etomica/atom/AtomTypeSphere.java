/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.simulation.Simulation;

/**
 * Atom type for a simple monatomic atom that has a length scale associated
 * with its size.  Position definition is the atom's coordinate 
 * (AtomPositionDefinitionSimple).
 */

public class AtomTypeSphere extends AtomTypeLeaf {
    
    double diameter, radius;
    
    public AtomTypeSphere(Simulation sim, AtomTypeGroup parentType) {
        this(parentType, sim.getDefaults().atomMass, sim.getDefaults().atomSize);
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
    public void setDiameter(double d) {
        if (d < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        diameter = d; radius = 0.5*d;
    }
}