/**
 * Defines units and dimensions that are used for conversions during input or output.  All quantities
 * used internally by the simulation are given in "simulation units", which is based on the
 * Angstrom as the unit of length, the Dalton as the unit of mass, and the picosecond as the unit
 * of time (a Dalton, or atomic mass unit (amu), is one gram divided by Avogadro's number).  The
 * classes given in this package can be used to convert between these units and another unit
 * that might be more natural or convenient.
 * The two key elements of the package are the Unit interface and the abstract Dimension class.
 * <h3>Dimension</h3>
 * Instances of Dimension can be used to indicate the physical nature of a value, e.g., whether
 * it is a quantity of mass, length, force, etc.  There are two places where this information is
 * routinely provided:
 * <ul>
 * <li>the package {@link etomica.data etomica.data} defines a Data class
 * which holds a DataInfo instance describing the Data; part of the information given in DataInfo
 * is a Dimension instance which specifies the physical dimensions of the data.
 * <li>the convention
 * of defining set/get methods to access and change fields of an object is supplemented with an annotation
 * on the get method specifying the dimension of the value.  For example, a class
 * that defines setSize and getSize methods to access the field <i>size</i>
 * (which represents, say, the diameter of a sphere) will also have an annotation on getSize indicating
 * dimensions of length.
 * </ul>
 * Dimensions are specified via a <i>signature</i>, which is an array of seven values indicating the
 * exponents of the fundamental dimensions that are combined to form the specified dimension.
 * The fundamental dimensions follow the convention defined by the
 * <a href="http://www.bipm.org/en/si/">SI system</a>, and are
 * <ul>
 * <li>length
 * <li>mass
 * <li>time
 * <li>electric current
 * <li>temperature
 * <li>number (of molecules)
 * <li>luminous intensity
 * </ul>
 * So, for example, the signature
 * of the energy dimension is length<sup>2</sup>-mass/time<sup>2</sup> and is given by the signature array (2, 1, -2, 0, 0, 0, 0).
 * The signature of energy/molecules (e.g., Joules/mole) is (2, 1, -2, 0, 0, -1, 0).
 * <p>
 * Dimension subclasses are defined for the fundamental dimensions and commonly encountered derived dimensions.
 * These classes have names such as Length, Time, Volume, Energy, and so on.  These classes all define
 * static singleton instances with the field name DIMENSION; the field SIM_UNIT in each class gives an
 * instance of Unit that corresponds to the unit derived from simulation units for that dimension.
 * Other dimensions can be defined by constructing instances of CompoundDimension.
 * <h3>Unit</h3>
 * Implementations of the Unit interface provide a convenient means to convert between
 * simulation units and some other particular unit, which need be done only when data is
 * read, written, or displayed.
 * <p>
 * SimpleUnit is a basic class that holds a conversion factor,
 * a Dimension instance and other descriptive information for implementing a basic Unit.  Many
 * specific units are defined by extending this class; examples include Kelvin, Bar, Joule, etc.
 * <p>
 * A PrefixedUnit class takes a Unit instance and a Prefix, which can be used to construct
 * units such as kilograms (combining Prefix.KILO with Gram.UNIT).
 * <p>
 * Derived units that are not already defined can be constructed using the CompoundUnit class.
 * <p>
 * The package {@link etomica.units.systems etomica.units.systems} defines constructs that
 * can collect the units defined by conventional unit systems, such as SI, cgs, English, etc.
  */
package etomica.units;