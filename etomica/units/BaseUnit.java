package etomica.units;

/**
 * Superclass for all base unit classes.  These classes provide a means for indicating
 * the physical units of a given quantity, and present methods for converting between units.
 * A BaseUnit can be used as is, or combined with a Prefix to form a
 * PrefixedUnit.  Units are employed by I/O classes (usually a Device or
 * Display) to handle unit conversions and labeling of graphic elements.
 * <br>
 * By convention, each subclass of any base unit will contain a static field named UNIT, which is a handle
 * to an instance of that unit.  One can access an instance of any unit class through this static member.
 * <br>
 * Each general base unit type (i.e., dimension) is defined as an abstract class.  Each of these abstract classes
 * contains an inner static subclass (named Sim) that defines the unit as derived from the basic simulation units (Dalton-A-ps)
 * Thus an instance of the base simulation units for any dimensioned quantity can be accessed by 
 * the handle BaseUnit.Energy.Sim.UNIT (e.g. for the energy unit).
 */

/* History
 * 03/11/04 (DAK) modified with changes introducing Unit interface, PrefixedUnit
 * class
 */
public abstract class BaseUnit implements Unit, java.io.Serializable {
    
    /**
     * Constructor defaulting to allow a prefix.
     * @param to conversion factor from this unit to simulation units
     * @param name string describing this unit
     * @param symbol symbol for this unit
     */
    public BaseUnit(double to, String name, String symbol, Dimension dimension) {
    	this(to, name, symbol, dimension, Prefix.ALLOWED);
    }
    public BaseUnit(double to, String name, String symbol, Dimension dimension, boolean prefixAllowed) {
    	setToSimConversionFactor(to);
    	this.name = name;
    	this.symbol = symbol;
    	this.dimension = dimension;
    	this.prefixAllowed = prefixAllowed;
    }
    
    /**
     * Conversion factor from the class unit to simulation units.
     * Set in constructor of subclass.
     */
    private double to;
    /**
     * Conversion factor from simulation units to the class unit
     * Set in constructor of subclass.
     */
    private double from;
    /**
     * A common name for the unit (e.g., Kelvins).  Written in plural.
     * Set in constructor of subclass.
     */
    private final String name;
    /**
     * A symbol for the unit (e.g., K for Kelvins)
     * Set in constructor of subclass.
     */
    private final String symbol;
    
    /**
     * Dimensions of this unit.
     */
    private final Dimension dimension;
    
    /**
     * Flag indicating whether setting a prefix (other than Null) is allowed.
     * Prefix is inappropriate, for example, with Angstrom unit, or a unit already defined
     * with a prefix, such as picosecond.  Default is true (prefix is allowed).
     * Value is modified appropriately in concrete subclass of Unit.
     */
    protected boolean prefixAllowed = true;       
    
    /**
     * Takes the given value in class units and converts it to simulation units.
     * @param x a value in units of this class
     * @return the value converted to simulation units
     */
    public final double toSim(double x) {return to*x;}
    
    /**
     * Takes the given value in simulation units and converts it to class units.
     * @param x a value in simulation units
     * @return the value converted to units of this class
     */
    public final double fromSim(double x) {return from*x;}
    
    /**
     * Allows subclass to change the conversion factor for this unit.  Used for
     * adjustable unit classes, such as Lennard-Jones units, for which the
     * conversion factor might be variable.
     * @param to  The new value of the conversion factor
     */
    protected void setToSimConversionFactor(double to) {
    	this.to = to;
    	from = 1.0/to;
    }
    
    /**
     * Accessor for common name of unit
     */
    public String toString() {return name;}
    
    /**
     * Accessor for symbol of unit
     */
    public String symbol() {return symbol;};
    
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
    public Dimension dimension() {return dimension;}
 
    //***** end of methods and fields of BaseUnit class *****//
    
    public interface D2 {  //marks a unit defined for a 2-dimensional space
    /**
     * FALSE_DEPTH is the size of a phony 3rd dimension ascribed to a 2-dimensional simulation
     * to permit its results to be reported in more familiar 3-dimensional units
     */
//        static final double FALSE_DEPTH = 5.0;  //Angstroms
    }  
    public interface D3 {}  //marks a unit defined for a 3-dimensional space (e.g., most pressure and volume units)

    /**
     * Null unit used for dimensionless quantities
     */
    public static class Null extends BaseUnit {
        public static final Null UNIT = new Null();
        public Null() {
        	super(1.0, "dimensionless", "", Dimension.NULL, Prefix.NOT_ALLOWED);
        }
    }

    /**
     * Undefined unit used for quantities with undefined or unknown units.
     */
    public static class Undefined extends BaseUnit {
        public static final Undefined UNIT = new Undefined();
        public Undefined() {
            super(1.0, "undefined", "", Dimension.UNDEFINED, Prefix.NOT_ALLOWED);
        }
    }

    /**
     * Simulation unit for the quantity of discrete things (e.g. molecules) is Count.
     * An example of another unit of this type is the mole.
     */
    public static abstract class Quantity extends BaseUnit {
		public Quantity(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.QUANTITY, prefixAllowed);}
		public Quantity(double to, String name, String symbol) {super(to, name, symbol, Dimension.QUANTITY);}
        public double to(Quantity u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Quantity u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Quantity {
            public static final Quantity UNIT = Count.UNIT;
			public Sim() {super(1.0, UNIT.toString(), UNIT.symbol(), UNIT.prefixAllowed());}
        }
    }
 
    /**
     * Simulation unit for a quantity representing a fractional amount
     *  (e.g. molecules) is Count.
     * An examplle of another unit of this type is the mole.
     */
    public static abstract class Fraction extends BaseUnit {
		public Fraction(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.FRACTION, prefixAllowed);}
		public Fraction(double to, String name, String symbol) {super(to, name, symbol, Dimension.FRACTION);}
        public double to(Quantity u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Quantity u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Fraction {
            public static final Fraction UNIT = Decimal.UNIT;
			public Sim() {super(1.0, UNIT.toString(), UNIT.symbol(), UNIT.prefixAllowed());}
        }
    }

    /**
     * Simulation unit of mass is the Dalton (1/AVOGADRO grams)
     */
    public static abstract class Mass extends BaseUnit {
		public Mass(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.MASS, prefixAllowed);}
		public Mass(double to, String name, String symbol) {super(to, name, symbol, Dimension.MASS);}
        public double to(Mass u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Mass u, double x) {return u.toSim(this.fromSim(x));}
                
        public static final class Sim extends Mass {
            public static final Mass UNIT = Dalton.UNIT;
			public Sim() {super(1.0, UNIT.toString(), UNIT.symbol(), UNIT.prefixAllowed());}
        }
    }

    /**
     * Simulation unit of length is the Angstrom
     */
    public static abstract class Length extends BaseUnit {
		public Length(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.LENGTH, prefixAllowed);}
		public Length(double to, String name, String symbol) {super(to, name, symbol, Dimension.LENGTH);}
        public double to(Length u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Length u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Length {
            public static final Length UNIT = Angstrom.UNIT;
            public static double TO_PIXELS = 10.0;
            public Sim() {super(1.0, UNIT.toString(), UNIT.symbol(), UNIT.prefixAllowed());}
        }
    }

    /**
     * Simulation unit of time is the picosecond
     */
    public static abstract class Time extends BaseUnit {
		public Time(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.TIME, prefixAllowed);}
		public Time(double to, String name, String symbol) {super(to, name, symbol, Dimension.TIME);}
        public double to(Time u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Time u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Time {
            public static final Time UNIT = Picosecond.UNIT;
			public Sim() {super(1.0, UNIT.toString(), UNIT.symbol(), UNIT.prefixAllowed());}
        }
    }
    
    /**
     * Simulation angular units are radians
     */
    public static abstract class Angle extends BaseUnit {
		public Angle(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.ANGLE, prefixAllowed);}
		public Angle(double to, String name, String symbol) {super(to, name, symbol, Dimension.ANGLE);}
        public double to(Angle u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Angle u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Angle {
            public static final Angle UNIT = Radian.UNIT;
			public Sim() {super(1.0, UNIT.toString(), UNIT.symbol(), UNIT.prefixAllowed());}
        }
    }
        
    /**
     * Simulation unit of charge is (D-A^3/ps^2)^(1/2)
     */
    public static abstract class Charge extends BaseUnit {
		public Charge(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.CHARGE, prefixAllowed);}
		public Charge(double to, String name, String symbol) {super(to, name, symbol, Dimension.CHARGE);}
        public double to(Charge u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Charge u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Charge {
            public static final Charge UNIT = new Sim();
			public Sim() {super(1.0, "sim charge units", "(D-A^3/ps^2)^(1/2)", Prefix.NOT_ALLOWED);}
        }
    }

    /**
     * Simulation unit of electrostatic dipole moment is (Daltons-A^5/ps^2)^(1/2)
     */
    public static abstract class Dipole extends BaseUnit {
		public Dipole(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.DIPOLE, prefixAllowed);}
		public Dipole(double to, String name, String symbol) {super(to, name, symbol, Dimension.DIPOLE);}
        public double to(Dipole u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Dipole u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Dipole {
            public static final Dipole UNIT = new Sim();
            public Sim() {super(1.0, "sim dipole units", "(D-A^5/ps^2)^(1/2)", Prefix.NOT_ALLOWED);}
		}
    }

    /**
     * Simulation unit of energy is D-A^2/ps^2
     */
    public static abstract class Energy extends BaseUnit {
		public Energy(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.ENERGY, prefixAllowed);}
		public Energy(double to, String name, String symbol) {super(to, name, symbol, Dimension.ENERGY);}
		private Energy(double to, String name, String symbol, Dimension dimension, boolean prefixAllowed) {super(to, name, symbol, dimension, prefixAllowed);}
		private Energy(double to, String name, String symbol, Dimension dimension) {super(to, name, symbol, dimension);}
        public double to(Energy u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Energy u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Energy {
            public static final Energy UNIT = new Sim();
            public Sim() {super(1.0, "sim energy units", "D-A^2/ps^2", Prefix.NOT_ALLOWED);}
        }
    }
    
    /**
     * Simulation unit of temperature is the simulation energy unit,
     * obtained by multiplying the temperature by Boltzmann's constant
     */
    public static abstract class Temperature extends Energy {
		public Temperature(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.TEMPERATURE, prefixAllowed);}
		public Temperature(double to, String name, String symbol) {super(to, name, symbol, Dimension.TEMPERATURE);}
        public static final class Sim extends Temperature {
            public static final Temperature UNIT = new Sim();
            public Sim() {super(1.0, "sim temperature units", "kB amu-A^2/ps^2", Prefix.NOT_ALLOWED);}
        }
    }
    
    /**
     * Simulation unit of (3D) pressure is (D-A/ps^2)/A^2 = D/(A-ps^2)
     */
    public static abstract class Pressure extends BaseUnit {
		public Pressure(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.PRESSURE, prefixAllowed);}
		public Pressure(double to, String name, String symbol) {super(to, name, symbol, Dimension.PRESSURE);}
		private Pressure(double to, String name, String symbol, Dimension dimension, boolean prefixAllowed) {super(to, name, symbol, dimension, prefixAllowed);}
		private Pressure(double to, String name, String symbol, Dimension dimension) {super(to, name, symbol, dimension);}
        public double to(Pressure u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Pressure u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Pressure {
            public static final Pressure UNIT = new Sim();
            public Sim() {super(1.0, "sim pressure units", "D/(A-ps^2)", Prefix.NOT_ALLOWED);}
        }
    }
    
    /**
     * Simulation unit of (2D) pressure is (amu-A/ps^2)/A = amu/ps^2
     */
    public static abstract class Pressure2D extends Pressure implements D2 {
		public Pressure2D(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.PRESSURE2D, prefixAllowed);}
		public Pressure2D(double to, String name, String symbol) {super(to, name, symbol, Dimension.PRESSURE2D);}
        public Dimension dimension() {return Dimension.PRESSURE2D;}
        //if converting to a "false" 3D pressure, divide by false depth
        public double to(Pressure u, double x) {
            return (u instanceof D3) ? u.fromSim(this.toSim(x/BaseUnitPseudo3D.FALSE_DEPTH)) : u.fromSim(this.toSim(x));
        }
        //if converting from a "false" 3D pressure, multiply by false depth
        public double from(Pressure u, double x) {
            return (u instanceof D3) ? u.toSim(this.fromSim(x*BaseUnitPseudo3D.FALSE_DEPTH)) : u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Pressure2D {
            public static final Pressure2D UNIT = new Sim();
            public Sim() {super(1.0, "sim 2-D pressure units", "D/ps^2", Prefix.NOT_ALLOWED);}
        }
    }
    
    /**
     * Simulation unit of (3D) volume is A^3
     */
    public static abstract class Volume extends BaseUnit {
		public Volume(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.VOLUME, prefixAllowed);}
		public Volume(double to, String name, String symbol) {super(to, name, symbol, Dimension.VOLUME);}
		public Volume(double to, String name, String symbol, Dimension dimension, boolean prefixAllowed) {super(to, name, symbol, dimension, prefixAllowed);}
		public Volume(double to, String name, String symbol, Dimension dimension) {super(to, name, symbol, dimension);}
        public double to(Volume u, double x) {return u.fromSim(this.toSim(x));}
        public double from(Volume u, double x) {return u.toSim(this.fromSim(x));}
        
        public static final class Sim extends Volume {
            public static final Volume UNIT = new Sim();
            public Sim() {super(1.0, "sim volume units", Angstrom.UNIT.symbol()+"^3", Prefix.NOT_ALLOWED);}
        }
    }
    
    /**
     * Simulation unit of (2D) volume is A^2
     */
    public static abstract class Volume2D extends Volume implements D2 {
		public Volume2D(double to, String name, String symbol, boolean prefixAllowed) {super(to, name, symbol, Dimension.VOLUME2D, prefixAllowed);}
		public Volume2D(double to, String name, String symbol) {super(to, name, symbol, Dimension.VOLUME2D);}
        //if converting to a "false" 3D volume, multiply by false depth
        public double to(Volume u, double x) {
            return (u instanceof D3) ? u.fromSim(this.toSim(x*BaseUnitPseudo3D.FALSE_DEPTH)) : u.fromSim(this.toSim(x));
        }
        //if converting from a "false" 3D volume, divide by false depth
        public double from(Volume u, double x) {
            return (u instanceof D3) ? u.toSim(this.fromSim(x/BaseUnitPseudo3D.FALSE_DEPTH)) : u.toSim(this.fromSim(x));}
       
        public static final class Sim extends Volume2D {
            public static final Volume2D UNIT = new Sim();
           public Sim() {super(1.0, "sim 2-D volume units", Angstrom.UNIT.symbol()+"^2", Prefix.NOT_ALLOWED);}
        }
    }
    
    //others to be defined include velocity, momentum, etc.
    
    /**
     * Returns an array of all available BaseUnit classes having the given dimension.
     * Finds them by performing instrospection of the classes in the units directory.
     */
    public static Class[] all(Dimension dimension) {
        if(dimension == null) throw new IllegalArgumentException("null argument for dimension passed to BaseUnit.all()");
        if(dimension == Dimension.NULL) return new Class[] {dimension.baseUnit()};
	    Class baseUnitClass = dimension.baseUnit();
	    java.io.File dir = new java.io.File(etomica.Default.CLASS_DIRECTORY+"/units");
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
	           String classname = "etomica.units."+files[i];
	           allClasses[i] = Class.forName(classname);
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	        if(allClasses[i]!=null && baseUnitClass.isAssignableFrom(allClasses[i])) someClasses[nClass++] = allClasses[i];
	    }
	    allClasses = new Class[nClass];
	    for(int i=0; i<nClass; i++) {allClasses[i] = someClasses[i];}
	    return allClasses;
    }
}
