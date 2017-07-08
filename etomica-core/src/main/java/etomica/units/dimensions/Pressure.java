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
 * Dimension for (3D) pressure. 
 */
public final class Pressure extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Pressure();
    /**
     * Simulation unit of pressure is (D-A/ps^2)/A^2 = D/(A-ps^2).
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION, 1, "sim pressure units", "D/(\u00c5-ps^2)", Prefix.NOT_ALLOWED);

    private Pressure() {
        super("Pressure", -1, 1, -2);
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        return unitSystem.pressure();
    }

    /**
     * Returns the Dimension instance for the "pressure" appropriate to the given
     * dimension. In particular, for D = 3, returns Pressure.DIMENSION; for D = 2,
     * returns Pressure2D.DIMENSION.
     * 
     * @throws IllegalArgumentException
     *             if D is not equal to 2, or 3.
     */
   public static Dimension dimension(int D) {
        switch(D) {
            case 3: 
                return Pressure.DIMENSION;
            case 2:
                return Pressure2D.DIMENSION;
//            case 1:
//                return Length.DIMENSION;
            default:
                throw new IllegalArgumentException("Pressure dimension defined only for D = 2, 3; you gave D = "+D);
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
