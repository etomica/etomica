/**
 * Provides classes that specify a system of units, such as SI, English, etc.  Unit classes are
 * sometimes used to convert values on input and output; the conversion is performed between
 * the units employed for internal calculations, and more recognizable units such at bars our Joules.
 * Typically a Unit will be determined from a Dimension instance associated with the I/O quantity.
 * The Dimension can provide an appropriate Unit when given a specific UnitSystem instance.
 * <p>
 * The packages {@link etomica.units etomica.units} and {@link etomica.units.dimensions etomica.units.dimensions} define the Unit and Dimension classes.
  */
package etomica.units.systems;