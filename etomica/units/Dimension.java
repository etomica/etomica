package etomica.units;
import etomica.Simulation;

/**
 * Collection of interfaces for dimensions of physical quantities.
 * These are used by get[FIELDNAME]Dimension() methods to determine the
 * physical dimensions of the field variable [FIELDNAME].  This information is used
 * by Devices and Displays to present options for setting the I/O units.<br>
 *
 * This interface defines a static final instance of each inner interface.
 * These fields are returned by the get[FIELDNAME]Dimension() methods to
 * specify the dimension of the quantity described by FIELDNAME.
 * For example, if the field timeStep has dimensions of time, then this would
 * be indicated by the method as follows:<br>
 * public Dimension getTimeStepDimension() {return Dimension.TIME;}
 * These instances are also returned by Meter and Modulator to indicate the dimensions of 
 * the quantity being measured or modulated.
 *
 * @author David Kofke
 */

//Java2 imports
//import java.util.LinkedList;

import etomica.utility.java2.LinkedList;

public abstract class Dimension implements java.io.Serializable {
    
    public Dimension() {}
    
    /**
     * Returns the unit corresponding to this dimension as given by the unitSystem object held by Simulation.
     * This becomes the default unit used for input/output in Devices and Displays.
     */
    public abstract Unit defaultIOUnit();
    /**
     * Class object for the base unit corresponding to this dimension.
     */
    public abstract Class baseUnit();
    
    /**
     * The signature is the exponents of each of the base dimensions forming the
     * given dimension.  Base dimensions are, in order: mass, length, time,
     * number.
     */
    public abstract double[] signature();
    
    /**
     * Dimension for a dimensionless quantity
     */
    public static class Null extends Dimension {
        private Null() {} //singleton; access via static instances, defined below
        static double[] signature = {0., 0., 0., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Dimensionless";}
        public Unit defaultIOUnit() {return Unit.NULL;}
        public Class baseUnit() {return BaseUnit.Null.class;}
    }
    
	/**
	 * Dimension indicating that a dimensions class is not yet defined
	 */
	public static class Undefined extends Dimension {
		private Undefined() {} //singleton; access via static instances, defined below
		static double[] signature = {0., 0., 0., 0.};
		public double[] signature() {return signature;}
		public String toString() {return "Undefined";}
		public Unit defaultIOUnit() {return Unit.NULL;}
		public Class baseUnit() {return BaseUnit.Null.class;}
	}

    /**
     * Dimension for a counted quantity, such as number of molecules
     */
    public static class Quantity extends Dimension {
        private Quantity() {}
        static double[] signature = {0., 0., 0., 1.};
        public double[] signature() {return signature;}
        public String toString() {return "Quantity";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().quantity();}
        public Class baseUnit() {return BaseUnit.Quantity.class;}
    }
    
    /**
     * Dimension for a fraction quantity, such as mole fraction.
     */
    public static class Fraction extends Dimension {
        private Fraction() {}
        static double[] signature = {0., 0., 0., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Decimal";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().fraction();}
        public Class baseUnit() {return BaseUnit.Fraction.class;}
   	
    }
    //remaining dimensions have obvious interpretations
    public static class Mass extends Dimension {
        private Mass() {}
        static double[] signature = {1., 0., 0., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Mass";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().mass();}
        public Class baseUnit() {return BaseUnit.Mass.class;}
    }
    public static class Length extends Dimension {
        private Length() {}
        static double[] signature = {0., 1., 0., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Length";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().length();}
        public Class baseUnit() {return BaseUnit.Length.class;}
    }
    public static class Time extends Dimension {
        private Time() {}
        static double[] signature = {0., 0., 1., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Time";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().time();}
        public Class baseUnit() {return BaseUnit.Time.class;}
    }
    public static class Angle extends Dimension {
        private Angle() {}
        static double[] signature = {0., 0., 0., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Angle";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().angle();}
        public Class baseUnit() {return BaseUnit.Angle.class;}
    }
    public static class Charge extends Dimension {//(D-A^3/ps^2)^(1/2)
        private Charge() {}
        static double[] signature = {0.5, 1.5, -1., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Charge";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().charge();}
        public Class baseUnit() {return BaseUnit.Charge.class;}
    }
    public static class Dipole extends Dimension {//(D-A^5/ps^2)^(1/2)
        private Dipole() {}
        static double[] signature = {0.5, 2.5, -1., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Dipole moment";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().dipole();}
        public Class baseUnit() {return BaseUnit.Dipole.class;}
    }
    public static class Energy extends Dimension {//D-A^2/ps^2
        private Energy() {}
        static double[] signature = {1., 2., -2., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Energy";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().energy();}
        public Class baseUnit() {return BaseUnit.Energy.class;}
    }
    public static class Temperature extends Energy {
        private Temperature() {}
        public String toString() {return "Temperature";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().temperature();}
        public Class baseUnit() {return BaseUnit.Temperature.class;}
    }
    public static class Pressure extends Dimension {//(D-A/ps^2)/A^2 = D/(A-ps^2)
        private Pressure() {}
        static double[] signature = {1., -1., -2., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Pressure";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().pressure(Simulation.instance.space().D());}
        public Class baseUnit() {return BaseUnit.Pressure.class;}
    }
    public static class Pressure2D extends Pressure {//(D-A/ps^2)/A = D/ps^2
        private Pressure2D() {}
        static double[] signature = {1., 0., -2., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "2D pressure";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().pressure(2);}
        public Class baseUnit() {return BaseUnit.Pressure2D.class;}
    }
    public static class Volume extends Dimension {
        private Volume() {}
        static double[] signature = {0., 3., 0., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "Volume";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().volume(Simulation.instance.space().D());}
        public Class baseUnit() {return BaseUnit.Volume.class;}
    }
    public static class Volume2D extends Volume {
        private Volume2D() {}
        static double[] signature = {0., 2., 0., 0.};
        public double[] signature() {return signature;}
        public String toString() {return "2D volume";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().volume(2);}
        public Class baseUnit() {return BaseUnit.Volume2D.class;}
    }
    
    //Singleton instances of each dimension
	public static final Dimension NULL = new Null();
    public static final Dimension QUANTITY = new Quantity();
	public static final Dimension FRACTION = new Fraction();
    public static final Dimension MASS = new Mass();
    public static final Dimension LENGTH = new Length();
    public static final Dimension TIME = new Time();
    public static final Dimension ANGLE = new Angle();
    public static final Dimension CHARGE = new Charge();
    public static final Dimension DIPOLE = new Dipole();
    public static final Dimension ENERGY = new Energy();
    public static final Dimension TEMPERATURE = new Temperature();
    public static final Dimension PRESSURE = new Pressure();
    public static final Dimension PRESSURE2D = new Pressure2D();
    public static final Dimension VOLUME = new Volume();
    public static final Dimension VOLUME2D = new Volume2D();
	public static final Dimension UNDEFINED = new Undefined();
   public static final Dimension[] ALL = new Dimension[] {
        NULL, QUANTITY, MASS, LENGTH, TIME, ANGLE, CHARGE, DIPOLE, ENERGY, 
        TEMPERATURE, PRESSURE, PRESSURE2D, VOLUME, VOLUME2D, UNDEFINED};

    
    /**
     * Returns all dimension classes with the same signature as the one given.
     */
    public static Dimension[] convertSignature(double[] sig) {
 //       java.util.ArrayList dimList = new java.util.ArrayList(5);
        LinkedList dimList = new LinkedList();
        for(int i=0; i<ALL.length; i++) {
            double[] dSig = ALL[i].signature();
            if(sig[0]==dSig[0] && sig[1]==dSig[1] && sig[2]==dSig[2] && sig[3]==dSig[3]) {
                dimList.add(ALL[i]);
            }
        }
        Dimension[] dimArr = new Dimension[dimList.size()];
        dimList.toArray(dimArr);
        return dimArr;
    }
    /**
     * Method to determine the dimension of a property via introspection.
     *
     * @param obj an instance of the object having the property of interest
     * @param property the name of the property, such that getPROPERTYDimension would return the dimension of the property
     * @param beaninfo BeanInfo class obtained by prior introspection of the object
     */
     //used by modulator, propertysheet
     public static Dimension introspect(Object obj, String property, java.beans.BeanInfo beaninfo) {
        java.beans.MethodDescriptor[] methods = beaninfo.getMethodDescriptors();
        for(int i=0; i<methods.length; i++) {
            String methodName = methods[i].getMethod().getName();
            if(methodName.equalsIgnoreCase("get"+property+"Dimension")) {
                try {
                    return (Dimension)methods[i].getMethod().invoke(obj, null);
                }
                catch(java.lang.reflect.InvocationTargetException ex) {
                    System.out.println("InvocationTargetException in Dimension.introspect");
                    return null;
                }
                catch(IllegalAccessException ex) {
                    System.out.println("IllegalAccessException in Dimension.introspect");
                    return null;
                }
                catch(ClassCastException ex) {
                    System.out.println("Method "+methodName+" does not return a Dimension object in Dimension.introspect");
                    return null;
                }
            }//end of if(methodName...) block
        }//end of for loop
        return NULL;
     }//end of introspect method
}