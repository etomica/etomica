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
 * Dimension for all volume units.
 */
public final class Volume extends Dimension {
    
    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Volume();
    
    /**
     * Simulation unit of volume is cubic Angstroms.
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "cubic Angstroms", "\u00c5^3", Prefix.NOT_ALLOWED);

    private Volume() {
        super("Volume", 3, 0, 0, 0, 0, 0, 0);
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.volume();
    }

    /**
     * Returns the Dimension instance for the "volume" appropriate to the given
     * dimension. In particular, for D = 3, returns Volume.DIMENSION; for D = 2,
     * returns Area.DIMENSION; for D = 1, returns Length.DIMENSION.
     * 
     * @throws IllegalArgumentException
     *             if D is not equal to 1, 2, or 3.
     */
   public static Dimension dimension(int D) {
        switch(D) {
            case 3: 
                return Volume.DIMENSION;
            case 2:
                return Area.DIMENSION;
            case 1:
                return Length.DIMENSION;
            default:
                throw new IllegalArgumentException("Volume dimension defined only for D = 1, 2, 3; you gave D = "+D);
        }
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
