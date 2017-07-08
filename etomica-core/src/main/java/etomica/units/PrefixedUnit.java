/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.units;

import etomica.units.dimensions.Dimension;

/**
 * Implementation of the Unit interface formed from a base Unit (e. g., grams) 
 * and a Prefix (e. g., kilo). The base unit cannot be changed
 * after construction. However, the prefix can be changed using the setPrefix
 * method.
 *
 * @see SimpleUnit
 * @see Prefix
 */

/* History
 * 03/11/04 (DAK) new, from Unit
 */

public class PrefixedUnit implements Unit, java.io.Serializable {
    
    /**
     * Scaling prefix for the baseUnit (e.g., milli, micro, kilo, etc.).
     * Conversion to the baseUnit takes into account this prefix.
     */
    private Prefix prefix;
    private final Unit baseUnit;
    private double to;
    private double from;
    private String name;
    private String symbol;
        
    /**
     * Constructs a PrefixedUnit using the given base and a Null prefix
     */
    public PrefixedUnit(Unit base) {
        this(Prefix.NULL, base);
    }
    /**
     * Constructs a PrefixedUnit using the given prefix and base
     */
    public PrefixedUnit(Prefix prefix, Unit base) {
        baseUnit = base;
        setPrefix(prefix);
    }
    
    /**
     * Returns the dimension of this baseUnit.
     * For example, the dimension of grams is mass.
     */
     public Dimension dimension() {return baseUnit.dimension();}
     
    /**
     * Changes prefix to given instance if prefixAllowed flag is true
     */
    public final void setPrefix(Prefix pre) {
        prefix = baseUnit.prefixAllowed() ? pre : Prefix.NULL;
        to = prefix.value()*baseUnit.toSim(1.0);
        from = 1.0/to;
        name = prefix.toString() + baseUnit.toString();
        symbol = prefix.symbol() + baseUnit.symbol();
    }
    
    public Prefix getPrefix() {return prefix;}
    
    /**
     * Returns the prefix of the baseUnit.  Equivalent to getPrefix().
     */
    public Prefix prefix() {return prefix;}
    /**
     * Returns the base baseUnit of the baseUnit.
     */
    public Unit unit() {return baseUnit;}
    
    /**
     * Takes the given value in class units (considering prefix) and converts it to simulation units.
     * @param x a value in units of this class
     * @return the value converted to simulation units
     */
    public final double toSim(double x) {return to*x;}
    
    /**
     * Takes the given value in simulation units and converts it to class units (considering prefix).
     * @param x a value in simulation units
     * @return the value converted to units of this class
     */
    public final double fromSim(double x) {return from*x;}
    
    /**
     * Accessor for common name of baseUnit
     */
    public String toString() {return name;}
    
    /**
     * Accessor for symbol of baseUnit
     */
    public String symbol() {return symbol;};
    
    /**
     * Returns false to indicate that a prefix cannot be applied to an already
     * prefixed baseUnit.
     * @see etomica.units.Unit#prefixAllowed()
     */
    public boolean prefixAllowed() {return false;}
    
    private static final long serialVersionUID = 1;


}
