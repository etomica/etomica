package etomica.units.systems;
import etomica.units.*;
/**
  * Meter-Kilogram-Second system of units
  */
    public class MKS extends UnitSystem {
    	
    	/**
    	 * Singleton instance of this unit system.
    	 */
    	public static final MKS SYSTEM = new MKS();
    	
    	private MKS() {}
    	
        public Unit quantity() {return Mole.UNIT;}
        public Unit fraction() {return Decimal.UNIT;}
        public Unit mass() {return new PrefixedUnit(Prefix.KILO, Gram.UNIT);}
        public Unit length() {return Meter.UNIT;}
        public Unit time() {return Second.UNIT;}
        public Unit angle() {return Radian.UNIT;}
        public Unit charge() {return Coulomb.UNIT;}
        public Unit dipole() {return Debye.UNIT;}  //??
        public Unit energy() {return Joule.UNIT;}
        public Unit temperature() {return Kelvin.UNIT;}
        public Unit pressure(int D) {return (D==2) ?
            (Unit)new BaseUnitPseudo3D.Pressure(Bar.UNIT) :
            (Unit)Bar.UNIT;}
        public Unit volume(int D) {return (D==2) ?
            (Unit)new BaseUnitPseudo3D.Volume(CubicMeter.UNIT) :
            (Unit)CubicMeter.UNIT;}
    }