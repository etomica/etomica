package simulate.units;

/**
 * Superclass for all base unit classes.  These classes provide a means for indicating
 * the physical units of a given quantity, and present methods for converting between units.
 * A BaseUnit is combined with a Prefix to form a Unit, which is then employed by I/O classes
 * (usually a Device or Display) to handle unit conversions and labeling of graphic elements.
 * <br>
 * By convention, each subclass of any base unit will contain a static field named UNIT, which is a handle
 * to an instance of that unit.  One can access an instance of any unit class through this static member.
 * <br>
 * Each general base unit type (i.e., dimension) is defined as an abstract class.  Each of these abstract classes
 * contains an inner static subclass (named Sim) that defines the unit as derived from the basic simulation units (Dalton-A-ps)
 * Thus an instance of the base simulation units for any dimensioned quantity can be accessed by 
 * the handle BaseUnit.Energy.Sim.UNIT (e.g. for the energy unit).  
 * <br>
 * Each base unit has a toPixels method that can be used to specify how a value of
 * a quantity in that unit would be scaled to render it some way on screen.
 * The scaling factor is declared static in the Sim subclass of that unit.
 */
public abstract class BaseUnit implements java.io.Serializable {
    
    public BaseUnit() {}  //constructor
    
    /**
     * Conversion factor from the class unit to simulation units.
     * Set in constructor of subclass.
     */
    protected double to = 1.0;
    /**
     * Conversion factor from simulation units to the class unit
     * Set in constructor of subclass.
     */
    protected double from = 1.0;
    /**
     * A common name for the unit (e.g., Kelvins).  Written in plural.
     * Set in constructor of subclass.
     */
    protected String name = "";
    /**
     * A symbol for the unit (e.g., K for Kelvins)
     * Set in constructor of subclass.
     */
    protected String symbol = "";
    
    /**
     * Flag indicating whether setting a prefix (other than Null) is allowed.
     * Prefix is inappropriate, for example, with Angstrom unit, or a unit already defined
     * with a prefix, such as picosecond.  Default is true (prefix is allowed).
     * Value is modified appropriately in concrete subclass of Unit.
     */
    protected boolean prefixAllowed = true;       
    
    /**
     * Takes the given value in class units (considering prefix) and converts it to simulation units.
     * @param x a value in units of this class
     * @return the value converted to simulation units
     */
    public final double toSim(double x) {return to*x;}
    
    /**
     * Takes the given value in simulation units and converts it to class units (considering prefix).
     * @param x a value in simulation units
     * @return the value converted to units of this class
     */
    public final double fromSim(double x) {return from*x;}
    
    /**
     * Accessor for common name of unit
     */
    public String toString() {return name;}
    
    /**
     * Accessor for symbol of unit
     */
    public String symbol() {return symbol;};

    /**
     * Convert given value from units of this class to display pixels.
     * Used perhaps if some representation of this value is displayed as a graphic to the screen.
     * Conversion factor for each dimension is held as a static constant in the Sim subclass 
     * of the general Unit subclass (e.g., Mass.Sim.TO_PIXELS is the factor for mass).
     */
    public abstract double toPixels(double x);
    
    /**
     * Returns flag indicating whether a prefix is allowed with this unit.
     * Some units (such as Angstroms) are not normally defined with a prefix attached,
     * and for such units this flag can be set to false to prohibit the application of a prefix.
     * This indication is usually made in the constructor of the base unit class.
     */
    public boolean prefixAllowed() {return prefixAllowed;}
     
    /**
     * Returns the dimension of this base unit.
     * For example, the dimension of grams is mass
     */
    public abstract Dimension dimension();
 
    //***** end of methods and fields of BaseUnit class *****//
    
    public interface D2 {  //marks a unit defined for a 2-dimensional space
    /**
     * FALSE_DEPTH is the size of a phony 3rd dimension ascribed to a 2-dimensional simulation
     * to permit its results to be reported in more familiar 3-dimensional units
     */
        static final double FALSE_DEPTH = 5.0;  //Angstroms
    }  
    public interface D3 {}  //marks a unit defined for a 3-dimensional space (e.g., most pressure and volume units)

    /**
     * Null unit used for dimensionless quantities
     */
    public static class Null extends BaseUnit {
        public Dimension dimension() {return Dimension.NULL;}
        public static final Null UNIT = new Null();
        public static double TO_PIXELS = 1.0;
        public Null() {
            prefixAllowed = false;
            name = "dimensionless";
        }
        public double toPixels(double x) {return x*TO_PIXELS;}
    }
    
    /**
     * Simulation unit for the quantity of discrete things (e.g. molecules) is Count.
     * An examplle of another unit of this type is the mole.
     */
    public static abstract class Quantity extends BaseUnit {
        public Dimension dimension() {return Dimension.QUANTITY;}
        public double to(Quantity u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Quantity u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        public static final class Sim extends Quantity {
            public static final Quantity UNIT = Count.UNIT;
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = UNIT.toString();  symbol = UNIT.symbol();}
        }
    }
    
    /**
     * Simulation unit of mass is the Dalton (1/AVOGADRO grams)
     */
    public static abstract class Mass extends BaseUnit {
        public Dimension dimension() {return Dimension.MASS;}
        public double to(Mass u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Mass u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Mass {
            public static final Mass UNIT = Dalton.UNIT;
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = UNIT.toString();  symbol = UNIT.symbol();}
        }
    }

    /**
     * Simulation unit of length is the Angstrom
     */
    public static abstract class Length extends BaseUnit {
        public Dimension dimension() {return Dimension.LENGTH;}
        public double to(Length u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Length u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Length {
            public static final Length UNIT = Angstrom.UNIT;
            public static double TO_PIXELS = 10.0;
            public Sim() {prefixAllowed = false; name = UNIT.toString();  symbol = UNIT.symbol();}
        }
    }

    /**
     * Simulation unit of time is the picosecond
     */
    public static abstract class Time extends BaseUnit {
        public Dimension dimension() {return Dimension.TIME;}
        public double to(Time u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Time u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Time {
            public static final Time UNIT = Picosecond.UNIT;
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = UNIT.toString();  symbol = UNIT.symbol();}
        }
    }
    
    /**
     * Simulation angular units are radians
     */
    public static abstract class Angle extends BaseUnit {
        public Dimension dimension() {return Dimension.ANGLE;}
        public double to(Angle u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Angle u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Angle {
            public static final Angle UNIT = Radian.UNIT;
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = true; name = UNIT.toString(); symbol = UNIT.symbol();}
        }
    }
        
    /**
     * Simulation unit of charge is (D-A^3/ps^2)^(1/2)
     */
    public static abstract class Charge extends BaseUnit {
        public Dimension dimension() {return Dimension.CHARGE;}
        public double to(Charge u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Charge u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Charge {
            public static final Charge UNIT = new Sim();
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = "sim charge units";  symbol = "(D-A^3/ps^2)^(1/2)";}
        }
    }

    /**
     * Simulation unit of electrostatic dipole moment is (Daltons-A^5/ps^2)^(1/2)
     */
    public static abstract class Dipole extends BaseUnit {
        public Dimension dimension() {return Dimension.DIPOLE;}
        public double to(Dipole u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Dipole u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Dipole {
            public static final Dipole UNIT = new Sim();
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = "sim dipole units";  symbol = "(D-A^5/ps^2)^(1/2)";}
        }
    }

    /**
     * Simulation unit of energy is D-A^2/ps^2
     */
    public static abstract class Energy extends BaseUnit {
        public Dimension dimension() {return Dimension.ENERGY;}
        public double to(Energy u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Energy u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Energy {
            public static final Energy UNIT = new Sim();
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = "sim energy units";  symbol = "D-A^2/ps^2";}
        }
    }
    
    /**
     * Simulation unit of temperature is the simulation energy unit,
     * obtained by multiplying the temperature by Boltzmann's constant
     */
    public static abstract class Temperature extends Energy {
        public Dimension dimension() {return Dimension.TEMPERATURE;}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        public static final class Sim extends Temperature {
            public static final Temperature UNIT = new Sim();
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = "times kB, sim units";  symbol = "kB amu-A^2/ps^2";}
        }
    }
    
    /**
     * Simulation unit of (3D) pressure is (D-A/ps^2)/A^2 = D/(A-ps^2)
     */
    public static abstract class Pressure extends BaseUnit {
        public Dimension dimension() {return Dimension.PRESSURE;}
        public double to(Pressure u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Pressure u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Pressure {
            public static final Pressure UNIT = new Sim();
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = "sim pressure units";  symbol = "D/(A-ps^2)";}
        }
    }
    
    /**
     * Simulation unit of (2D) pressure is (amu-A/ps^2)/A = amu/ps^2
     */
    public static abstract class Pressure2D extends Pressure implements D2 {
        public Dimension dimension() {return Dimension.PRESSURE2D;}
        //if converting to a "false" 3D pressure, divide by false depth
        public double to(Pressure u, double x) {
            return (u instanceof D3) ? u.fromSim(this.toSim(x/BaseUnit.D2.FALSE_DEPTH)) : u.fromSim(this.toSim(x));
        }
        //if converting from a "false" 3D pressure, multiply by false depth
        public double from(Pressure u, double x) {
            return (u instanceof D3) ? u.toSim(this.fromSim(x*BaseUnit.D2.FALSE_DEPTH)) : u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Pressure2D {
            public static final Pressure2D UNIT = new Sim();
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = "sim 2-D pressure units";  symbol = "D/ps^2";}
        }
    }
    
    /**
     * Simulation unit of (3D) volume is A^3
     */
    public static abstract class Volume extends BaseUnit {
        public Dimension dimension() {return Dimension.VOLUME;}
        public double to(Volume u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Volume u, double x) {return u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Volume {
            public static final Volume UNIT = new Sim();
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = "sim volume units";  symbol = Angstrom.UNIT.symbol()+"^3";}
        }
    }
    
    /**
     * Simulation unit of (2D) volume is A^2
     */
    public static abstract class Volume2D extends Volume implements D2 {
        public Dimension dimension() {return Dimension.VOLUME2D;}
        //if converting to a "false" 3D volume, multiply by false depth
        public double to(Volume u, double x) {
            return (u instanceof D3) ? u.fromSim(this.toSim(x*BaseUnit.D2.FALSE_DEPTH)) : u.fromSim(this.toSim(x));
        }
        //if converting from a "false" 3D volume, divide by false depth
        public double from(Volume u, double x) {
            return (u instanceof D3) ? u.toSim(this.fromSim(x/BaseUnit.D2.FALSE_DEPTH)) : u.toSim(this.fromSim(x));}
        public double toPixels(double x) {return this.toSim(x)*Sim.TO_PIXELS;}
        
        public static final class Sim extends Volume2D {
            public static final Volume2D UNIT = new Sim();
            public static double TO_PIXELS = 1.0;
            public Sim() {prefixAllowed = false; name = "sim 2-D volume units";  symbol = Angstrom.UNIT.symbol()+"^2";}
        }
    }
    
    //others to be defined include velocity, momentum, etc.
    
    /**
     * Returns an array of all available BaseUnit classes having the given dimension.
     * Finds them by performing instrospection of the classes in the units directory.
     */
    public static Class[] all(Dimension dimension) {
        if(dimension == null) throw new IllegalArgumentException("null argument for dimension passed to BaseUnit.all()");
	    Class baseUnitClass = dimension.baseUnit();
	    java.io.File dir = new java.io.File(simulate.Default.CLASS_DIRECTORY+"/units");
	    String[] files = dir.list(new java.io.FilenameFilter() {
	        public boolean accept(java.io.File d, String name) {
	              return !name.startsWith("BaseUnit") ||
	                        name.endsWith("Sim.class");}
	        });
//	    java.util.Arrays.sort(files);  //won't autojar with this included
	    Class[] allClasses = new Class[files.length];
	    Class[] someClasses = new Class[files.length];
	    int nClass = 0;
	    for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");  //drop the ".class" suffix
	        files[i] = files[i].substring(0,idx);
	        allClasses[i] = null;
	        try{
	           String classname = "simulate.units."+files[i];
	           allClasses[i] = Class.forName(classname);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	        if(allClasses[i]!=null && baseUnitClass.isAssignableFrom(allClasses[i])) someClasses[nClass++] = allClasses[i];
	    }
	    allClasses = new Class[nClass];
	    for(int i=0; i<nClass; i++) {allClasses[i] = someClasses[i];}
	    return allClasses;
    }
}