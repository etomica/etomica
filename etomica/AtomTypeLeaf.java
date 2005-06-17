/*
 * Created on May 20, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica;

import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomPositionDefinitionSimple;
import etomica.units.Dimension;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class AtomTypeLeaf extends AtomType {

    /**
     * @param parentType
     * @param positionDefinition
     */
    public AtomTypeLeaf(AtomType parentType, double mass) {
        this(parentType, new AtomPositionDefinitionSimple(), mass);
     }
    /**
     * @param parentType
     * @param positionDefinition
     * @param mass
     */
    public AtomTypeLeaf(AtomType parentType,
            AtomPositionDefinition positionDefinition, double mass) {
        super(parentType, positionDefinition);
        setMass(mass);
    }
    
    /**
     * Returns true, indicating that this is a leaf type.
     */
    public boolean isLeaf() {
        return true;
    }
    
    /**
     * Sets  mass of this atom and updates reciprocal mass accordingly.  Setting
     * mass to largest machine double (Double.MAX_VALUE) causes reciprocal mass 
     * to be set to zero.
     * 
     * @param mass   new value for mass
     */
    public void setMass(double m) {
        mass = m;
        rm = (m==Double.MAX_VALUE) ? 0.0 : 1.0/mass;
    }
    public final double getMass() {return mass;}
    public final double rm() {return rm;}
    public final Dimension getMassDimension() {return Dimension.MASS;}
    
    public double mass, rm;

}
