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
 * Base for all area units. Simulation unit of area is A^2
 */
public class Area extends Dimension {
    
    public static final Dimension DIMENSION = new Area();
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "square Angstroms", "\u00c5^2", Prefix.NOT_ALLOWED);

    private Area() {
        super("Area", 2, 0, 0, 0, 0, 0, 0);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.area();
//        return new CompoundUnit(new Unit[] {unitSystem.length()}, new double[] {2.0});
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
