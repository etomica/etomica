/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.systems;

import etomica.units.Angstrom;
import etomica.units.dimensions.*;
import etomica.units.Candela;
import etomica.units.Count;
import etomica.units.Dalton;
import etomica.units.Decimal;
import etomica.units.Picosecond;
import etomica.units.Radian;
import etomica.units.Unit;

/**
 * Specifies a system of units.  Defines methods that return specific units for
 * various types of quantities.
 */
public abstract class UnitSystem  {

    /**
     * @return unit defining an amount in this unit system, typically a direct count or a mole
     */
    public abstract Unit quantity();

    /**
     * @return unit defining a fractional amount in this unit system, typically a decimal or percent
     */
    public abstract Unit fraction();

    /**
     * @return mass unit defined in this unit system
     */
    public abstract Unit mass();

    /**
     * @return length unit defined in this unit system
     */
    public abstract Unit length();

    /**
     * @return time unit defined in this unit system
     */
    public abstract Unit time();

    /**
     * @return angle unit defined in this unit system, typically radians or degrees
     */
    public abstract Unit angle();

    /**
     * @return electric-charge unit defined in this unit system
     */
    public abstract Unit charge();

    /**
     * @return electric-current unit defined in this unit system
     */
    public abstract Unit current();

    /**
     * @return electric-dipole unit defined in this unit system
     */
    public abstract Unit dipole();

    /**
     * @return force unit defined in this unit system
     */
    public abstract Unit force();

    /**
     * @return energy unit defined in this unit system
     */
    public abstract Unit energy();

    /**
     * @return power (energy/time) unit defined in this unit system
     */
    public abstract Unit power();

    /**
     * @return temperature unit defined in this unit system
     */
    public abstract Unit temperature();

    /**
     * @return pressure unit defined in this unit system
     */
    public abstract Unit pressure();

    /**
     * @return volume unit defined in this unit system
     */
    public abstract Unit volume();

    /**
     * @return area unit defined in this unit system
     */
    public abstract Unit area();

    /**
     * @return viscosity unit defined in this unit system
     */
    public abstract Unit viscosity();

    /**
     * @return luminous intensity unit defined in this unit system
     */
    public Unit luminousIntensity() {
        return Candela.UNIT;
    }

    /**
     * Base dimensions are {length, mass, time, current, temperature, quantity, luminous intensity}
     * @return a list of the base units defined for this unit system
     */
    public Unit[] baseUnits() {
        return new Unit[] {length(), mass(), time(), current(), temperature(), quantity(), luminousIntensity()};
    }
    
 /**
  * System of units based on simulation units of Daltons, picoseconds, and Angstroms
  */
    public static class Sim extends UnitSystem {
        
        public Unit quantity() {return Count.UNIT;}
        public Unit fraction() {return Decimal.UNIT;}
        public Unit mass() {return Dalton.UNIT;}
        public Unit length() {return Angstrom.UNIT;}
        public Unit time() {return Picosecond.UNIT;}
        public Unit angle() {return Radian.UNIT;}
        public Unit charge() {return Charge.SIM_UNIT;}
        public Unit current() {return Current.SIM_UNIT;}
        public Unit dipole() {return Dipole.SIM_UNIT;}
        public Unit force() {return Force.SIM_UNIT;}
        public Unit energy() {return Energy.SIM_UNIT;}
        public Unit power() {return Power.SIM_UNIT;}
        public Unit temperature() {return Temperature.SIM_UNIT;}
        public Unit pressure() {return Pressure.SIM_UNIT;}
        public Unit pressure2D() {return Pressure2D.SIM_UNIT;}
        public Unit volume() {return Volume.SIM_UNIT;}
        public Unit area() {return Area.SIM_UNIT;}
        public Unit viscosity() {return Viscosity.SIM_UNIT;}
    }

    /**
     * Instance of simulation units, the unit system used for internal calculations.
     */
    public static final UnitSystem SIM = new Sim();
}
