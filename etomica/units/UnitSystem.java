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
        
        public Unit quantity() {return Count.UNIT;}
        public Unit mass() {return Dalton.UNIT;}
        public Unit length() {return Angstrom.UNIT;}
        public Unit time() {return Picosecond.UNIT;}
        public Unit angle() {return Radian.UNIT;}
        public Unit charge() {return BaseUnit.Charge.Sim.UNIT;}
        public Unit dipole() {return BaseUnit.Dipole.Sim.UNIT;}
        public Unit energy() {return BaseUnit.Energy.Sim.UNIT;}
        public Unit temperature() {return BaseUnit.Temperature.Sim.UNIT;}
        public Unit pressure(int D) {return (D==2) ? BaseUnit.Pressure2D.Sim.UNIT : BaseUnit.Pressure.Sim.UNIT;}
        public Unit volume(int D) {return (D==2) ? BaseUnit.Volume2D.Sim.UNIT : BaseUnit.Volume.Sim.UNIT;}
    }
    public static final UnitSystem SIM = new Sim();
    //will also define CGS, atomic, English
}