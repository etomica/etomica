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
 * Dimension for all units of viscosity, Mass/(Length-Time)
 */
public final class Viscosity extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Viscosity();
    /**
     * The simulation unit is Dalton/(A-ps).
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "daltons/(angstrom-ps)", "D/(\u00c5-ps)", Prefix.NOT_ALLOWED);

    private Viscosity() {
        super("Viscosity", -1, 1, -1);// LMTCtNl;
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.mass();
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
