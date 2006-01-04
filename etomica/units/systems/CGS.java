package etomica.units.systems;

import java.io.ObjectStreamException;

import etomica.units.Barye;
import etomica.units.CompoundUnit;
import etomica.units.Coulomb;
import etomica.units.CubicCentimeter;
import etomica.units.Decimal;
import etomica.units.Dyne;
import etomica.units.Erg;
import etomica.units.Gram;
import etomica.units.Kelvin;
import etomica.units.Meter;
import etomica.units.Mole;
import etomica.units.Poise;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Radian;
import etomica.units.Second;
import etomica.units.Unit;

/**
 * Centimeter-Gram-Second system of units.
 */
public class CGS extends UnitSystem {

    /**
     * Singleton instance of this unit system.
     */
    public static final CGS SYSTEM = new CGS();

    private CGS() {
    }

    public Unit quantity() {
        return Mole.UNIT;
    }

    public Unit fraction() {
        return Decimal.UNIT;
    }

    public Unit mass() {
        return Gram.UNIT;
    }

    public Unit length() {
        return CENTIMETER;
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
        return COULOMB_CENTIMETER;
    }

    public Unit force() {
        return Dyne.UNIT;
    }
    
    public Unit energy() {
        return Erg.UNIT;
    }

    public Unit temperature() {
        return Kelvin.UNIT;
    }

    public Unit pressure() {
        return Barye.UNIT;
    }

    public Unit volume() {
        return CubicCentimeter.UNIT;
    }
    
    public Unit area() {
        return SQUARE_CENTIMETER;
    }
    
    public Unit viscosity() {
        return Poise.UNIT;
    }
    
    private static final Unit CENTIMETER = new PrefixedUnit(Prefix.CENTI, Meter.UNIT);
    private static final Unit SQUARE_CENTIMETER = new CompoundUnit(new Unit[] {CENTIMETER}, new double[] {2});
    private static final Unit COULOMB_CENTIMETER = 
        new CompoundUnit(new Unit[] {Coulomb.UNIT, CENTIMETER}, new double[] {1, 1});
    
    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton SYSTEM
     */
    private Object readResolve() throws ObjectStreamException {
        return SYSTEM;
    }
    
    private static final long serialVersionUID = 1;

}