/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.api.IAtomTypeSphere;
import etomica.api.IElement;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;

/**
 * Atom type for a simple monatomic atom that has a length scale associated
 * with its size.  Position definition is the atom's coordinate 
 * (AtomPositionDefinitionSimple).
 */

public class AtomTypeSphere extends AtomTypeLeaf implements IAtomTypeSphere {
    
    private static final long serialVersionUID = 1L;
    protected double diameter;
    
    public AtomTypeSphere(Simulation sim) {
        this(new ElementSimple(sim), 1.0);
    }
    
    public AtomTypeSphere(IElement element) {
        this(element, 1.0);
    }
    
    public AtomTypeSphere(IElement element, double d) {
        super(element);
        setDiameter(d);
    }
                
    /* (non-Javadoc)
	 * @see etomica.atom.IAtomTypeSphere#getDiameter()
	 */
    public double getDiameter() {return diameter;}
    
    /* (non-Javadoc)
	 * @see etomica.atom.IAtomTypeSphere#setDiameter(double)
	 */
    public void setDiameter(double d) {
        if (d < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        diameter = d;
    }
}