package etomica.units;

/**
 * Base class for classes describing the default units to be used for I/O operations.
 * An instance of a UnitSystem class is held in a static field in the Simulation class 
 * to define the default system of units.
 *
 * @see Unit
 * @see Dimension
 * @see etomica.Simulation
 */
public abstract class UnitSystem implements java.io.Serializable {
    
    public abstract Unit quantity();
    public abstract Unit mass();
    public abstract Unit length();
    public abstract Unit time();
    public abstract Unit angle();
    public abstract Unit charge();
    public abstract Unit dipole();
    public abstract Unit energy();
    public abstract Unit temperature();
    public abstract Unit pressure(int D);
    public abstract Unit volume(int D);
    
    //Default pressure and volume units are those for 2-D system
    public Unit pressure() {return pressure(2);}
    public Unit volume() {return volume(2);}

 /**
  * System of units based on simulation units of Daltons, picoseconds, and Angstroms
  */
    public static class Sim extends UnitSystem {
        
        public Unit quantity() {return new Unit(Count.UNIT);}
        public Unit mass() {return new Unit(Dalton.UNIT);}
        public Unit length() {return new Unit(Angstrom.UNIT);}
        public Unit time() {return new Unit(Picosecond.UNIT);}
        public Unit angle() {return new Unit(Radian.UNIT);}
        public Unit charge() {return new Unit(BaseUnit.Charge.Sim.UNIT);}
        public Unit dipole() {return new Unit(BaseUnit.Dipole.Sim.UNIT);}
        public Unit energy() {return new Unit(BaseUnit.Energy.Sim.UNIT);}
        public Unit temperature() {return new Unit(BaseUnit.Temperature.Sim.UNIT);}
        public Unit pressure(int D) {return (D==2) ? new Unit(BaseUnit.Pressure2D.Sim.UNIT) : new Unit(BaseUnit.Pressure.Sim.UNIT);}
        public Unit volume(int D) {return (D==2) ? new Unit(BaseUnit.Volume2D.Sim.UNIT) : new Unit(BaseUnit.Volume.Sim.UNIT);}
    }
 /**
  * Meter-Kilogram-Second system of units
  */
    public static class MKS extends UnitSystem {
        public Unit quantity() {return new Unit(Mole.UNIT);}
        public Unit mass() {return new Unit(Prefix.KILO, Gram.UNIT);}
        public Unit length() {return new Unit(Meter.UNIT);}
        public Unit time() {return new Unit(Second.UNIT);}
        public Unit angle() {return new Unit(Radian.UNIT);}
        public Unit charge() {return new Unit(Coulomb.UNIT);}  //need to define Coulomb
        public Unit dipole() {return new Unit(Debye.UNIT);}  //??
        public Unit energy() {return new Unit(Joule.UNIT);}
        public Unit temperature() {return new Unit(Kelvin.UNIT);}
        public Unit pressure(int D) {return new Unit(Bar.UNIT);}  //??
        public Unit volume(int D) {return new Unit(CubicMeter.UNIT);} //??
    }
    
    /**
     * Lennard-Jones units.
     * All quantities are made dimensionless with respect to a characteristic
     * size (sigma), energy (epsilon) and mass.
     * Values of sigma, epsilon, and mass are set via the corresponding static
     * set methods in etomica.units.LennardJones class.
     */
     public static class LJ extends UnitSystem {        
        public Unit quantity() {return new Unit(Count.UNIT);}
        public Unit mass() {return new Unit(LennardJones.Mass.UNIT);}
        public Unit length() {return new Unit(LennardJones.Length.UNIT);}
        public Unit time() {return new Unit(LennardJones.Time.UNIT);}
        public Unit angle() {return new Unit(Radian.UNIT);}
        public Unit charge() {return new Unit(LennardJones.Charge.UNIT);}  
        public Unit dipole() {return new Unit(LennardJones.Dipole.UNIT);}  
        public Unit energy() {return new Unit(LennardJones.Energy.UNIT);}
        public Unit temperature() {return new Unit(LennardJones.Energy.UNIT);}
        public Unit pressure(int D) {return (D==2) ? new Unit(LennardJones.Pressure2D.UNIT) : new Unit(LennardJones.Pressure.UNIT);}  
        public Unit volume(int D) {return (D==2) ? new Unit(LennardJones.Volume2D.UNIT) : new Unit(LennardJones.Volume.UNIT);} 
        
     }   
   
    public static final UnitSystem MKS = new MKS();
    public static final UnitSystem LJ = new LJ();
    public static final UnitSystem SIM = new Sim();
    //will also define CGS, atomic, English
}