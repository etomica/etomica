/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;

/**
 * Unit formed from other units raised to specified powers. Used to construct unit when desired unit is not
 * specifically implemented in an existing class. For example, a surface-tension unit in N/m could
 * be created by constructing this class as
 * <pre>{@code
 * new CompoundUnit(new Unit[] {Newton.UNIT, Meter.UNIT}, new double[] {1.0,-1.0});
 * }</pre>
 */
public class CompoundUnit implements Unit {

    private final double from;
    private final String symbol;
    private final Dimension dimension;

    /**
     * Constructs unit with default name and abbreviation constructed from given units and exponents
     *
     * @param units     array of units used to specify this unit
     * @param exponents array exponents associated with each unit.
     * @throws IllegalArgumentException if argument arrays are of different length
     */
    public CompoundUnit(Unit[] units, double[] exponents) {
        this(units, exponents, makeName(units, exponents), makeSymbol(units, exponents));
    }

    /**
     * @param units     list of units used to specify this unit
     * @param exponents exponents associated with each unit
     * @param name      a common name used to reference this unit (e.g., Newton)
     * @param symbol    an abbreviated name used to reference this unit (e.g., N)
     * @throws IllegalArgumentException if units and exponents arrays are of different length
     */
    public CompoundUnit(Unit[] units, double[] exponents, String name, String symbol) {
        if (units.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundUnit constructor must be arrays of the same length.  Given were arrays of length " + units.length + " and " + exponents.length);
        }
        this.symbol = symbol;
        double f = 1.0;
        Dimension[] dimensions = new Dimension[units.length];
        for (int i = 0; i < units.length; i++) {
            f *= Math.pow(units[i].fromSim(1.0), exponents[i]);
            dimensions[i] = units[i].dimension();
        }
        from = f;
        dimension = new CompoundDimension(dimensions, exponents);
    }

    // used by constructor
    private static String makeName(Unit[] units, double[] exponents) {
        if (units.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundUnit constructor must be arrays of the same length.  Given were arrays of length " + units.length + " and " + exponents.length);
        }
        String name = "";
        for (int i = 0; i < units.length; i++) {
            if (exponents[i] == 0.0) continue;
            name += "(" + units[i].toString();
            if (exponents[i] != 1.0) {
                name += "^" + exponents[i];
            }
            name += ")";
        }
        return name;
    }

    // used by constructor
    private static String makeSymbol(Unit[] units, double[] exponents) {
        if (units.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundUnit constructor must be arrays of the same length.  Given were arrays of length " + units.length + " and " + exponents.length);
        }
        String symbol = "";
        for (int i = 0; i < units.length; i++) {
            if (exponents[i] == 0.0) continue;
            symbol += units[i].symbol();
            if (exponents[i] != 1.0) {
                symbol += "^" + exponents[i];
            }
            if (i != units.length - 1) symbol += "-";
        }
        return symbol;
    }

    public Dimension dimension() {
        return dimension;
    }

    public double toSim(double x) {
        return x / from;
    }

    public double fromSim(double x) {
        return x * from;
    }

    public String symbol() {
        return symbol;
    }

    public boolean prefixAllowed() {
        return false;
    }
}

