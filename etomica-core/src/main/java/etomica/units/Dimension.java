/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;
import java.util.Arrays;

import etomica.units.systems.UnitSystem;

/**
 * Parent of all Dimension classes, which describe the physical dimensions
 * (e.g., mass, length, force) of a quantity.  All dimensions are derived from
 * seven base dimensions (length, mass, time, electrical current, temperature,
 * number or mole, and luminosity).  Dimension is coded as a "signature", which
 * is an array of seven value that are the exponents of these base dimensions
 * forming the specified dimension.
 */
public class Dimension implements java.io.Serializable {

    /**
     * Number of base dimensions, equal to seven.  
     * Base dimensions are: length, mass, time, current, temperature, number, luminosity.
     */
    public static final int N_BASE = 7;

    public Dimension(String name, double length, double mass, double time) {
        this(name, length, mass, time, 0, 0, 0, 0);
    }
    
    public Dimension(String name, double length, double mass, double time, double current,
            double temperature, double number, double luminosity) {
        this(name, new double[] {length, mass, time, current, temperature, number, luminosity});
    }
    
    public Dimension(String name, double[] signature) {
        if(signature.length != N_BASE) {
            throw new IllegalArgumentException("Incorrect length of signature array given to Dimension constructor. Given value = "+signature.length+"; expected value: "+N_BASE);
        }
        this.signature = signature.clone();
        this.name = name;
    }
    
    /**
     * Returns the unit of this dimension as derived in the given system of units.
     * Default constructs a CompoundUnit using the base units of the unit system.
     */
    public Unit getUnit(UnitSystem unitSystem) {
        return new CompoundUnit(unitSystem.baseUnits(), signature());
    }
    
    public String toString() {
        return name;
    }
    /**
     * The signature is the exponents of each of the base dimensions forming the
     * given dimension.  Base dimensions are, in order: length, mass, time, current, 
     * temperature, number, luminosity.
     */
    public double[] signature() {
        return signature.clone();
    }
        
    /**
     * Returns true if the given object is a Dimension instance with the same signature as this.
     */
    public boolean equals(Object dim) {
        if(dim instanceof Dimension) {
            return Arrays.equals(this.signature, ((Dimension)dim).signature);
        }
        return false;
    }
    
    private final double[] signature;
    private final String name;
    private static final long serialVersionUID = 1L;

    /**
     * Dimension used to indicate that a group of values are not all of
     * the same dimension.  Singleton.  The equals method for this instance
     * returns true only if its argument is the same instance.
     */
    public static Dimension MIXED = new Mixed();

//   public static final Dimension[] ALL = new Dimension[] {
//        NULL, QUANTITY, MASS, LENGTH, TIME, ANGLE, CHARGE, DIPOLE, ENERGY, 
//        TEMPERATURE, PRESSURE, PRESSURE2D, VOLUME, VOLUME2D, UNDEFINED};
//
    
//    /**
//     * Returns all dimension classes with the same signature as the one given.
//     */
//    public static Dimension[] convertSignature(double[] sig) {
// //       java.util.ArrayList dimList = new java.util.ArrayList(5);
//        LinkedList dimList = new LinkedList();
//        for(int i=0; i<ALL.length; i++) {
//            double[] dSig = ALL[i].signature();
//            if(sig[0]==dSig[0] && sig[1]==dSig[1] && sig[2]==dSig[2] && sig[3]==dSig[3]) {
//                dimList.add(ALL[i]);
//            }
//        }
//        Dimension[] dimArr = new Dimension[dimList.size()];
//        dimList.toArray(dimArr);
//        return dimArr;
//    }
    /**
     * Method to determine the dimension of a property via introspection.
     *
     * @param obj an instance of the object having the property of interest
     * @param property the name of the property, such that getPROPERTYDimension would return the dimension of the property
     * @param beaninfo BeanInfo class obtained by prior introspection of the object
     */
     //used by modifier, propertysheet
     public static Dimension introspect(Object obj, String property, java.beans.BeanInfo beaninfo) {
        java.beans.MethodDescriptor[] methods = beaninfo.getMethodDescriptors();
        for(int i=0; i<methods.length; i++) {
            String methodName = methods[i].getMethod().getName();
            if(methodName.equalsIgnoreCase("get"+property+"Dimension")) {
                try {
                    return (Dimension)methods[i].getMethod().invoke(obj);
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
        return Null.DIMENSION;
     }//end of introspect method
     
     private static class Mixed extends Dimension {
         
         private Mixed() {
             super("Mixed", 
                 Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
         }
         /**
          * Required to guarantee singleton when deserializing.
          * 
          * @return the singleton MIXED
          */
         private Object readResolve() {
             return MIXED;
         }
         
         public boolean equals(Object object) {
             return (this == object);
         }
         
         public Unit getUnit(UnitSystem unitSystem) {
             return Undefined.UNIT;
         }
         
         private static final long serialVersionUID = 1;
     }
}
