/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Dimension;

/**
 * Interface used to specify the physical units to be used when inputting or
 * outputting a quantity. All internal calculations are performed using
 * simulation units, which are all derived from the picosecond, Dalton, and
 * Angstrom.  Classes implementing this interface convert the simulation units
 * to a more common unit for interfacing with the user.
 * <br>
 * I/O classes (normally Device and Display) have associated with them a
 * class that implements the Unit interface.  In addition to providing
 * conversions, the unit class also provides textual labels that can be used by
 * a Device or Display to indicate the units is it using.
 */
public interface Unit {

	/**
	 * Returns the dimension of the unit. 
	 * For example, the dimension of grams is mass.
	 */
	 public Dimension dimension();

	/**
	 * Takes the given value in class units (considering prefix) and converts it to simulation units.
	 * @param x a value in units of the class
	 * @return the value converted to simulation units
	 */
	public double toSim(double x);
    
	/**
	 * Takes the given value in simulation units and converts it to class units (considering prefix).
	 * @param x a value in simulation units
	 * @return the value converted to units of the class
	 */
	public double fromSim(double x);
    
	/**
	 * Returns the symbol of unit, such as "g" for grams.
	 */
	public String symbol();
	
	/**
	 * Indicates if the unit is appropriate for use with a prefix.
	 */
	public boolean prefixAllowed();

}
