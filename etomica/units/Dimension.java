package simulate.units;
import simulate.Simulation;

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
     * Dimension for a dimensionless quantity
     */
    public static class Null extends Dimension {
        public String toString() {return "Dimensionless";}
        public Unit defaultIOUnit() {return new Unit(BaseUnit.Null.UNIT);}
        public Class baseUnit() {return BaseUnit.Null.class;}
    }
    /**
     * Dimension for a counted quantity, such as number of molecules
     */
    public static class Quantity extends Dimension {
        public String toString() {return "Quantity";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().quantity();}
        public Class baseUnit() {return BaseUnit.Quantity.class;}
    }
    
    //remaining dimensions have obvious interpretations
    public static class Mass extends Dimension {
        public String toString() {return "Mass";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().mass();}
        public Class baseUnit() {return BaseUnit.Mass.class;}
    }
    public static class Length extends Dimension {
        public String toString() {return "Length";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().length();}
        public Class baseUnit() {return BaseUnit.Length.class;}
    }
    public static class Time extends Dimension {
        public String toString() {return "Time";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().time();}
        public Class baseUnit() {return BaseUnit.Time.class;}
    }
    public static class Angle extends Dimension {
        public String toString() {return "Angle";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().angle();}
        public Class baseUnit() {return BaseUnit.Angle.class;}
    }
    public static class Charge extends Dimension {
        public String toString() {return "Charge";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().charge();}
        public Class baseUnit() {return BaseUnit.Charge.class;}
    }
    public static class Dipole extends Dimension {
        public String toString() {return "Dipole moment";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().dipole();}
        public Class baseUnit() {return BaseUnit.Dipole.class;}
    }
    public static class Energy extends Dimension {
        public String toString() {return "Energy";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().energy();}
        public Class baseUnit() {return BaseUnit.Energy.class;}
    }
    public static class Temperature extends Energy {
        public String toString() {return "Temperature";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().temperature();}
        public Class baseUnit() {return BaseUnit.Temperature.class;}
    }
    public static class Pressure extends Dimension {
        public String toString() {return "Pressure";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().pressure(Simulation.instance.space().D());}
        public Class baseUnit() {return BaseUnit.Pressure.class;}
    }
    public static class Pressure2D extends Pressure {
        public String toString() {return "2D pressure";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().pressure(2);}
        public Class baseUnit() {return BaseUnit.Pressure2D.class;}
    }
    public static class Volume extends Dimension {
        public String toString() {return "Volume";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().volume(Simulation.instance.space().D());}
        public Class baseUnit() {return BaseUnit.Volume.class;}
    }
    public static class Volume2D extends Volume {
        public String toString() {return "2D volume";}
        public Unit defaultIOUnit() {return Simulation.unitSystem().volume(2);}
        public Class baseUnit() {return BaseUnit.Volume2D.class;}
    }
    
    //Convenience instances of each dimension
    public static final Dimension NULL = new Null();
    public static final Dimension QUANTITY = new Quantity();
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
        return null;
     }
}