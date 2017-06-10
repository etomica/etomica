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
 * The dimension for electrical current. Internally, charge q always represents
 * q/(4 pi epsilon0)^(1/2), so that the force via Coulomb's law is given by q1q2/r^2. 
 * This gives a quantity for charge having units of sqrt(force)*length, or mass^(1/2)-length^(3/2)/time, 
 * so the simulation unit of for the quantity representing charge is sqrt(D-A^3)/ps.
 * Consequently the simulation unit for current (charge/time) is sqrt(D-A^3)/ps^2.
 */
public class Current extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Current();
    /**
     * Simulation unit of electrical current.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1.0, "sim current units", "(D-A^3)^(1/2)/ps/(4 pi \u03B50)^(1/2)", Prefix.NOT_ALLOWED);

    private Current() {
        super("Current", 0, 0, 0, 1, 0, 0, 0);
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
