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
 * The dimension for electrical charge. Internally, charge q always represents
 * q/(4 pi epsilon0)^(1/2), so that the force via Coulomb's law is given by q1*q2/r^2. 
 * Internally, charge is a quantity having units of sqrt(force)*length, or mass^(1/2)-length^(3/2)/time, 
 * so the simulation unit for the quantity representing charge is sqrt(D-A^3)/ps.
 */
public class Charge extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Charge();
    /**
     * Simulation unit of charge.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1.0, "sim charge units", "(D-\u00c5^3/ps^2)^(1/2)(4\u03C0\u03B50)^(1/2)", Prefix.NOT_ALLOWED);//unicode is Angstroms, and 4 pi epsilon0

    private Charge() {
        super("Charge", 0, 0, 1, 1, 0, 0, 0);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.charge();
    }
}
