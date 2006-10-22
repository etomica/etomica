/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;

/**
 * Atom type for a simple monatomic atom that has a length scale associated
 * with its size.  Position definition is the atom's coordinate 
 * (AtomPositionDefinitionSimple).
 */

public class AtomTypeSphere extends AtomTypeLeaf {
    
    double diameter;
    
    public AtomTypeSphere(Simulation sim) {
        this(new ElementSimple(sim), sim.getDefaults().atomSize);
    }
    
    public AtomTypeSphere(Simulation sim, Element element) {
        this(element, sim.getDefaults().atomSize);
    }
    
    public AtomTypeSphere(Element element, double d) {
        super(element, new AtomPositionDefinitionSimple());
        setDiameter(d);
    }
                
    public double getDiameter() {return diameter;}
    
    /**
    * Sets diameter of this atom and updates radius accordingly.
    *
    * @param d   new value for diameter
    */
    public void setDiameter(double d) {
        if (d < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        diameter = d;
    }
}