/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.ISimulation;

/**
 * Atom type for a simple monatomic atom that has a length scale associated
 * with its size.  Position definition is the atom's coordinate 
 * (AtomPositionDefinitionSimple).
 */

public class AtomTypeSphere extends AtomTypeLeaf {
    
    private static final long serialVersionUID = 1L;
    protected double diameter;
    
    public AtomTypeSphere(ISimulation sim) {
        this(new ElementSimple(sim), 1.0);
    }
    
    public AtomTypeSphere(Element element) {
        this(element, 1.0);
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