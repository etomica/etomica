/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units.dimensions;

import etomica.units.CompoundUnit;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;

public class CompoundDimension extends Dimension {

    public CompoundDimension(Dimension[] dimensions, double[] exponents) {
        this(dimensions, exponents, makeName(dimensions, exponents));
    }
    
    public CompoundDimension(Dimension[] dimensions, double[] exponents, String name) {
        super(name, makeSignature(dimensions, exponents));
        this.dimensions = dimensions.clone();
        this.exponents = exponents.clone();
    }
    
    // used by constructor
    private static double[] makeSignature(Dimension[] dimensions, double[] exponents) {
        if(dimensions.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundDimension constructor must be arrays of the same length.  Given were arrays of length "+dimensions.length+" and "+exponents.length);
        }
        double[] sig = new double[N_BASE];
        for(int i=0; i<dimensions.length; i++) {
            for(int j=0; j<sig.length; j++) {
                sig[j] += dimensions[i].signature()[j]*exponents[i];
            }
        }
        return sig;
    }
    // used by constructor
    private static String makeName(Dimension[] dimensions, double[] exponents) {
        if(dimensions.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundDimension constructor must be arrays of the same length.  Given were arrays of length "+dimensions.length+" and "+exponents.length);
        }
        String symbol = "";
        for(int i=0; i<dimensions.length; i++) {
            symbol += "("+dimensions[i].toString()+"^"+exponents[i]+")";
        }
        return symbol;
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        Unit[] units = new Unit[dimensions.length];
        for(int i=0; i<units.length; i++) {
            units[i] = dimensions[i].getUnit(unitSystem);
        }
        return new CompoundUnit(units, exponents);
    }
    
    private final Dimension[] dimensions;
    private final double[] exponents;
    private static final long serialVersionUID = 1;
}
