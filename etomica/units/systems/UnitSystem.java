package etomica.units.systems;

import etomica.units.Angstrom;
import etomica.units.Area;
import etomica.units.Candela;
import etomica.units.Charge;
import etomica.units.Count;
import etomica.units.Dalton;
import etomica.units.Decimal;
import etomica.units.Dipole;
import etomica.units.Energy;
import etomica.units.Picosecond;
import etomica.units.Pressure;
import etomica.units.Pressure2D;
import etomica.units.Radian;
import etomica.units.Temperature;
import etomica.units.Unit;
import etomica.units.Volume;

/**
 * Specifies a system of units.  Defines methods that return specific units for
 * various types of quanities.
 */
public abstract class UnitSystem implements java.io.Serializable {
    
    public abstract Unit quantity();
    public abstract Unit fraction();
    public abstract Unit mass();
    public abstract Unit length();
    public abstract Unit time();
    public abstract Unit angle();
    public abstract Unit charge();
    public abstract Unit dipole();
    public abstract Unit energy();
    public abstract Unit temperature();
    public abstract Unit pressure();
    public abstract Unit volume();
    public abstract Unit area();
    public Unit luminousIntensity() {
        return Candela.UNIT;
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
        public Unit dipole() {return Dipole.SIM_UNIT;}
        public Unit energy() {return Energy.SIM_UNIT;}
        public Unit temperature() {return Temperature.SIM_UNIT;}
        public Unit pressure() {return Pressure.SIM_UNIT;}
        public Unit pressure2D() {return Pressure2D.SIM_UNIT;}
        public Unit volume() {return Volume.SIM_UNIT;}
        public Unit area() {return Area.SIM_UNIT;}
        private static final long serialVersionUID = 1;
    }
    public static final UnitSystem SIM = new Sim();
}