/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.systems;

import etomica.units.*;

/**
 * Meter-Kilogram-Second system of units
 */
public class MKS extends UnitSystem {

    public MKS() {
    }

    public Unit quantity() {
        return Mole.UNIT;
    }

    public Unit fraction() {
        return Decimal.UNIT;
    }

    public Unit mass() {
        return KILOGRAM;
    }

    public Unit length() {
        return Meter.UNIT;
    }

    public Unit time() {
        return Second.UNIT;
    }

    public Unit angle() {
        return Radian.UNIT;
    }

    public Unit charge() {
        return Coulomb.UNIT;
    }
    
    public Unit current() {
        return Ampere.UNIT;
    }

    public Unit dipole() {
        return COULOMB_METER;
    }

    public Unit force() {
        return Newton.UNIT;
    }
    
    public Unit energy() {
        return Joule.UNIT;
    }

    public Unit power(){return Watt.UNIT;}

    public Unit temperature() {
        return Kelvin.UNIT;
    }

    public Unit pressure() {
        return Bar.UNIT;
    }

    public Unit volume() {
        return CubicMeter.UNIT;
    }
    
    public Unit area() {
        return SQUARE_METER;
    }
    
    public Unit viscosity() {
        return PASCAL_SECOND;
    }
    
    private static final Unit SQUARE_METER = new CompoundUnit(new Unit[] {Meter.UNIT}, new double[] {2});
    private static final Unit KILOGRAM = new PrefixedUnit(Prefix.KILO, Gram.UNIT);
    private static final Unit COULOMB_METER = 
        new CompoundUnit(new Unit[] {Coulomb.UNIT, Meter.UNIT}, new double[] {1, 1});
    private static final Unit PASCAL_SECOND = 
        new CompoundUnit(new Unit[] {Pascal.UNIT, Second.UNIT}, new double[] {1, 1});
    
    private static final long serialVersionUID = 1;

}