/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import etomica.units.*;
import etomica.units.systems.UnitSystem;

/**
 * Dimension for all units of spatial frequency, 1/Length. Number of waves per unit length.
 */
public final class Wavenumber extends Dimension {

    /**
     * Singleton instance of this class.
     */
    public static final Dimension DIMENSION = new Wavenumber();
    /**
     * The simulation unit is 1/A
     */
    public static final Unit SIM_UNIT = new SimpleUnit(DIMENSION,1.0,"inverse Angstroms","\u00c5^-1",Prefix.NOT_ALLOWED);


    private Wavenumber() {
        super("Wavenumber", -1, 0, 0);// LMTCtNl;
    }

    public Unit getUnit(UnitSystem unitSystem) {
        return new CompoundUnit(new Unit[] {unitSystem.length()},new double[] {-1, 0,0});
    }

}
