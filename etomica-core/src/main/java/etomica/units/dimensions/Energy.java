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
 * Dimension for all energy units. Simulation unit of energy is D-A^2/ps^2
 */
public final class Energy extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Energy();
    /**
     * The simulation unit of energy is D-A^2/ps^2.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1.0, "sim energy units", "D-A^2/ps^2", Prefix.NOT_ALLOWED);

    private Energy() {
        super("Energy", 2, 1, -2);
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.energy();
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
