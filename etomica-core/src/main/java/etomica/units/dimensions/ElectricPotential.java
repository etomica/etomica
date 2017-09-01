/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import java.io.ObjectStreamException;

import etomica.units.Prefix;
import etomica.units.SimpleUnit;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * The dimension for electrostatic potential. The dimensions of electric
 * potential is energy/charge, which has fundamental dimensions of (length)^2
 * (mass) (time)^-3 (current)^-1.
 * <p>
 * Internally, charge is represented as q/(4 pi epsilon0)^1/2, and is a quantity
 * having units of mass^(1/2)-length^(3/2)/time, so the simulation unit for the
 * quantity representing electric potential is [(A)^2 (D) (ps)^-2] [(D)^-1/2
 * A^-3/2 ps] = A^(1/2) D^(1/2) ps^(-1).
 * 
 */
public class ElectricPotential extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new ElectricPotential();
    /**
     * Simulation unit of charge.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1.0, "sim electric-potential units", "(D-\u00c5)^(1/2)/ps/(4\03C0\u03B50)^(1/2)", Prefix.NOT_ALLOWED);//unicode is Angstroms, and 4 pi epsilon0

    private ElectricPotential() {
        super("Electric Potential", 2, 1, -3, -1, 0, 0, 0);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.charge();
    }

    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton DIMENSION
     */
    private Object readResolve() throws ObjectStreamException {
        return DIMENSION;
    }

    private static final long serialVersionUID = 1;

}
