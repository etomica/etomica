/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Jul 25, 2004 by kofke
 */
package etomica.math;

/**
 * @author kofke
 *
 * Data construct that can be used to define a double range.  Holds a
 * minimum and maximum value and can test if a given double is within
 * the range.  Each instance is immutable.
 */
public final class DoubleRange implements java.io.Serializable {

	private final double min;
	private final double max;
	/**
	 * Constructs the range defined by given values.  Minimum of range
	 * will be given by smaller of two values (if they differ) and
	 * maximum by larger.
	 */
	public DoubleRange(double bound1, double bound2) {
		if(bound1 < bound2) {
			min = bound1;
			max = bound2;
		} else {
			min = bound2;
			max = bound1;
		}
	}
	
	/**
	 * The minimium defined by this range.
	 * @return the value of the minimum
	 */
	public double minimum() {return min;}
	
	/**
	 * The maximum defined by this range.
	 * @return the value of the maximum.
	 */
	public double maximum() {return max;}
	
	/**
	 * Indicates if a value falls within the range of this instance. 
	 * @param k the tested value
	 * @return true if the given value between the range's 
	 * minimum and maximum values, inclusive.
	 */
	public boolean contains(double k) {
		return (k >= min) && (k <=max);
	}
	
	/**
	 * Tests for equality of mininum and maximum values of this 
	 * range with given range
	 * @param range
	 * @return true if minimum and maximum match.
	 */
	public boolean equals(DoubleRange range) {
		return (min==range.min && max==range.max);
	}
	
	/**
	 * Returns String of the form [min,max]
	 */
	public String toString() {
		return ("[" + min + "," + max + "]");
	}
}
