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
    
    public abstract Unit quantity();
    public abstract Unit fraction();
    public abstract Unit mass();
    public abstract Unit length();
    public abstract Unit time();
    public abstract Unit angle();
    public abstract Unit charge();
    public abstract Unit current();
    public abstract Unit dipole();
    public abstract Unit force();
    public abstract Unit energy();
    public abstract Unit power();
    public abstract Unit temperature();
    public abstract Unit pressure();
    public abstract Unit volume();
    public abstract Unit area();
    public abstract Unit viscosity();
    public Unit luminousIntensity() {
        return Candela.UNIT;
    }
    
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
    public static final UnitSystem SIM = new Sim();
}
