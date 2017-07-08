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
 * The temperature dimension. Internally, temperature always represents kT,
 * that is, the temperature multiplied by Boltzmann's constant; this gives a
 * quantity having dimensions of energy, and thus is in units of D-A^2/ps^2.  However,
 * temperature is treated as fundamental when defining its Dimension.
 */
public final class Temperature extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Temperature();
    /**
     * The simulation unit for temperature, which as kT has dimensions
     * of energy and is D-A^2/ps^2.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1.0, "sim temperature units", "kB D-A^2/ps^2", Prefix.NOT_ALLOWED);

    private Temperature() {
        super("Temperature", 0, 0, 0, 0, 1, 0, 0);// LMTCtNl
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.temperature();
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
