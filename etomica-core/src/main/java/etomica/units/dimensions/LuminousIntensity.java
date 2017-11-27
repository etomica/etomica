/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import java.io.ObjectStreamException;

import etomica.units.Candela;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

/**
 * The dimension for luminous intensity.  Simulation unit is the Candela.
 */
public final class LuminousIntensity extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new LuminousIntensity();
    /**
     * The Candela unit.
     */
    public static final Unit SIM_UNIT = Candela.UNIT;

    private LuminousIntensity() {
        super("Luminous Intensity", 0, 0, 0, 0, 0, 0, 1);// LMTCtNl
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.luminousIntensity();
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
