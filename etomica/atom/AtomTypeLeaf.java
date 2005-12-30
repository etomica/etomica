package etomica.atom;

import etomica.units.Dimension;
import etomica.units.Mass;

/**
 * Type for an atom that is a leaf in the species hierarchy. An atom of this
 * type is typically representing a physical atom, rather than a group of other
 * atoms.
 * 
 * @author andrew
 */

/*
 * Created on May 20, 2005
 */

public class AtomTypeLeaf extends AtomType {

    /**
     * Constructs type with position defined by AtomPositionDefinitionSimple.
     */
    public AtomTypeLeaf(AtomTypeGroup parentType, double mass) {
        this(parentType, new AtomPositionDefinitionSimple(), mass);
    }

    /**
     * Invokes parent constructor with and sets the mass to the given value.
     */
    public AtomTypeLeaf(AtomTypeGroup parentType,
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
     * Sets mass of this atom and updates reciprocal mass accordingly. Setting
     * mass to largest machine double (Double.MAX_VALUE) causes reciprocal mass
     * to be set to zero.
     * 
     * @param mass
     *            new value for mass
     */
    public void setMass(double m) {
        mass = m;
        rm = (m == Double.MAX_VALUE) ? 0.0 : 1.0 / mass;
    }

    /**
     * Returns the value of the mass.
     */
    public final double getMass() {
        return mass;
    }

    /**
     * Returns the reciprocal of the mass, 1.0/mass
     */
    public final double rm() {
        return rm;
    }

    /**
     * Returns Dimension.MASS, indicating that "mass" has dimensions of mass.
     */
    public final Dimension getMassDimension() {
        return Mass.DIMENSION;
    }

    public double mass, rm;

}