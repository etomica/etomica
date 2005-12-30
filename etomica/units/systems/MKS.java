package etomica.units.systems;

import java.io.ObjectStreamException;

import etomica.units.Bar;
import etomica.units.CompoundUnit;
import etomica.units.Coulomb;
import etomica.units.CubicMeter;
import etomica.units.Debye;
import etomica.units.Decimal;
import etomica.units.Gram;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Meter;
import etomica.units.Mole;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Radian;
import etomica.units.Second;
import etomica.units.Unit;
import etomica.units.UnitSystem;

/**
 * Meter-Kilogram-Second system of units
 */
public class MKS extends UnitSystem {

    /**
     * Singleton instance of this unit system.
     */
    public static final MKS SYSTEM = new MKS();

    private MKS() {
    }

    public Unit quantity() {
        return Mole.UNIT;
    }

    public Unit fraction() {
        return Decimal.UNIT;
    }

    public Unit mass() {
        return new PrefixedUnit(Prefix.KILO, Gram.UNIT);
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

    public Unit dipole() {
        return Debye.UNIT;
    } // ??

    public Unit energy() {
        return Joule.UNIT;
    }

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
    
    private static final Unit SQUARE_METER = new CompoundUnit(new Unit[] {Meter.UNIT}, new double[] {2});
    
    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton UNIT
     */
    private Object readResolve() throws ObjectStreamException {
        return SYSTEM;
    }
    
    private static final long serialVersionUID = 1;

}