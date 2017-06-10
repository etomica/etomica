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
 * Base unit for electrical dipole moment. Simulation unit is electron-Angstrom.
 */
public class Dipole extends Dimension {

    public static final Dimension DIMENSION = new Dipole();
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "sim dipole units", "e-\u00c5", Prefix.NOT_ALLOWED);

    private Dipole() {
        super("Dipole", 1, 0, 1, 1, 0, 0, 0);
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.dipole();
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
