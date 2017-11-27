/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Dimension;

/**
 * Superclass for all base unit classes. These classes provide a means for
 * indicating the physical units of a given quantity, and present methods for
 * converting between units. A BaseUnit can be used as is, or combined with a
 * Prefix to form a PrefixedUnit. Units are employed by I/O classes (usually a
 * Device or Display) to handle unit conversions and labeling of graphic
 * elements. <br>
 * By convention, each subclass of any base unit will contain a static field
 * named UNIT, which is a handle to an instance of that unit. One can access an
 * instance of any unit class through this static member. <br>
 * Each general base unit type (i.e., dimension) is defined as an abstract
 * class. Each of these abstract classes contains an inner static subclass
 * (named Sim) that defines the unit as derived from the basic simulation units
 * (Dalton-A-ps) Thus an instance of the base simulation units for any
 * dimensioned quantity can be accessed by the handle BaseUnit.Energy.Sim.UNIT
 * (e.g. for the energy unit).
 */
public class SimpleUnit implements Unit {

    /**
     * @param toSim
     *            conversion factor from this unit to simulation units
     * @param name
     *            string describing this unit
     * @param symbol
     *            symbol for this unit
     * @param prefixAllowed
     *            flag indicating if this unit is suitable for use with a prefix
     *            (e.g., kilo, nano)
     */
    public SimpleUnit(Dimension dimension, double toSim, String name, String symbol, boolean prefixAllowed) {
        this.dimension = dimension;
        fromSim = 1.0 / toSim;
        this.name = name;
        this.symbol = symbol;
        this.prefixAllowed = prefixAllowed;
    }

    /**
     * Returns the dimension of this base unit. For example, the dimension of
     * grams is mass.
     */
    public Dimension dimension() {
        return dimension;
    }

    /**
     * Takes the given value in class units and converts it to simulation units.
     * 
     * @param x
     *            a value in units of this class
     * @return the value converted to simulation units
     */
    public final double toSim(double x) {
        return x / fromSim;
    }

    /**
     * Takes the given value in simulation units and converts it to class units.
     * 
     * @param x
     *            a value in simulation units
     * @return the value converted to units of this class
     */
    public final double fromSim(double x) {
        return x * fromSim;
    }

    /**
     * @return the common name of this unit
     */
    public String toString() {
        return name;
    }

    /**
     * @return a symbol of this unit
     */
    public String symbol() {
        return symbol;
    };

    /**
     * Some units (such as Angstroms) are not normally defined with a prefix
     * attached, and for such units this flag can be set to false to prohibit
     * the application of a prefix. This indication is usually made in the
     * constructor of the base unit class.
     *
     * @return flag indicating whether a prefix is allowed with this unit
     */
    public boolean prefixAllowed() {
        return prefixAllowed;
    }

    /**
     * Conversion factor from simulation units to the class unit. Set in
     * constructor of subclass.
     */
    private double fromSim;

    /**
     * A common name for the unit (e.g., Kelvins). Written in plural. Set in
     * constructor of subclass.
     */
    private final String name;

    /**
     * A symbol for the unit (e.g., K for Kelvins) Set in constructor of
     * subclass.
     */
    private final String symbol;

    /**
     * Flag indicating whether setting a prefix (other than Null) is allowed.
     * Prefix is inappropriate, for example, with Angstrom unit, or a unit
     * already defined with a prefix, such as picosecond. Default is true
     * (prefix is allowed). Value is modified appropriately in concrete subclass
     * of Unit.
     */
    private final boolean prefixAllowed;
    
    /**
     * The physical dimension of the unit, e.g., mass, length, time.
     */
    private final Dimension dimension;

    // ***** end of methods and fields of SimpleUnit class *****//

//    /**
//     * Returns an array of all available BaseUnit classes having the given
//     * dimension. Finds them by performing instrospection of the classes in the
//     * units directory.
//     */
//    public static Class[] all(Dimension dimension) {
//        if (dimension == null)
//            throw new IllegalArgumentException("null argument for dimension passed to BaseUnit.all()");
//        if (dimension == Null.DIMENSION)
//            return new Class[] { dimension.getBaseUnitClass() };
//        Class baseUnitClass = dimension.getBaseUnitClass();
//        java.io.File dir = null;// new
//        // java.io.File(etomica.Default.CLASS_DIRECTORY+"/units");
//        String[] files = dir.list(new java.io.FilenameFilter() {
//            public boolean accept(java.io.File d, String name) {
//                return !name.startsWith("BaseUnit") || name.endsWith("Sim.class");
//            }
//        });
//        // java.util.Arrays.sort(files); //won't autojar with this included
//        Class[] allClasses = new Class[files.length];
//        Class[] someClasses = new Class[files.length];
//        int nClass = 0;
//        for (int i = 0; i < files.length; i++) {
//            int idx = files[i].lastIndexOf("."); // drop the ".class" suffix
//            files[i] = files[i].substring(0, idx);
//            allClasses[i] = null;
//            try {
//                String classname = "etomica.units." + files[i];
//                allClasses[i] = Class.forName(classname);
//            } catch (ClassNotFoundException e) {
//                System.out.println("Failed for " + files[i]);
//            }
//            if (allClasses[i] != null && baseUnitClass.isAssignableFrom(allClasses[i]))
//                someClasses[nClass++] = allClasses[i];
//        }
//        allClasses = new Class[nClass];
//        for (int i = 0; i < nClass; i++) {
//            allClasses[i] = someClasses[i];
//        }
//        return allClasses;
//    }
}
