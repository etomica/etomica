package simulate.units;

/**
 * Class used to specify the physical units to be used when inputting or outputting a quantity.
 * A Unit is usually formed from a BaseUnit (e.g., grams) and a Prefix (e.g., kilo), which together
 * determine the factor that is needed to convert from this unit to/from the simulation units. 
 * (All internal calculations are performed using simulation units, which are all derived from 
 *  the picosecond, Dalton, and Angstrom).
 * <br>
 * I/O classes (normally Device and Display) have associated with them a Unit class.  In addition 
 * to providing conversions, the unit class also provides textual labels that can be used by a 
 * Device or Display to indicate the units is it using.
 *
 * Once the unit is formed, the base unit cannot be changed.  However, the prefix can be changed
 * using the setPrefix method.
 *
 * @see UnitSystem
 * @see UnitRatio
 * @see BaseUnit
 * @see Prefix
 */
public class Unit implements java.io.Serializable {
    
    /**
     * Scaling prefix for the unit (e.g., milli, micro, kilo, etc.).
     * Conversion to the unit takes into account this prefix.
     */
    private Prefix prefix;
    private BaseUnit baseUnit;
    protected double to;
    protected double from;
    protected double baseTo;
    protected String name;
    protected String symbol;
    protected String baseName;
    protected String baseSymbol;
    protected Dimension dimension;
    private boolean prefixAllowed;
    
    /**
     * Constructs a Unit using the given base and a Null prefix
     */
    public Unit(BaseUnit base) {
        this(Prefix.NULL, base);
    }
    /**
     * Constructs a Unit using the given prefix and base
     */
    public Unit(Prefix p, BaseUnit base) {
        this(p, base.toSim(1.0), base.toString(), base.symbol(), base.prefixAllowed(), base.dimension());
        baseUnit = base;
    }
    /**
     * Constructs a Unit without a base unit, instead with explicit specification of all information usually provided by the base
     * 
     * @param pre the prefix for this unit
     * @param to  the factor for converting from this unit to simulation units
     * @param name the full name for this unit (e.g., Joule)
     * @param symbol an abbreviation for this unit (e.g., J)
     * @param prefixAllowed a flag indicating if this unit may take a non-Null prefix
     * @param dim the dimension of this unit (e.g., energy)
     */
    public Unit(Prefix pre, double to, String name, String symbol, boolean prefixAllowed, Dimension dim) {
        prefix = prefix;
        baseTo = to;
        baseName = name;
        baseSymbol = symbol;
        this.prefixAllowed = prefixAllowed;
        dimension = dim;
        setPrefix(pre);
    }
    
    /**
     * Returns the dimension of this unit.
     * For example, the dimension of grams is mass.
     */
     public Dimension dimension() {return dimension;}
     
    /**
     * Changes prefix to given instance if prefixAllowed flag is true
     */
    public final void setPrefix(Prefix pre) {
        prefix = prefixAllowed ? pre : Prefix.NULL;
        to = prefix.value()*baseTo;
        from = 1.0/to;
        name = prefix.toString() + baseName;
        symbol = prefix.symbol() + baseSymbol;
    }
    
    public Prefix getPrefix() {return prefix;}
    
    /**
     * Returns the prefix of the unit.  Equivalent to getPrefix().
     */
    public Prefix prefix() {return prefix;}
    /**
     * Returns the base unit of the unit.
     */
    public BaseUnit baseUnit() {return baseUnit;}
    
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
     * Accessor for common name of unit
     */
    public String toString() {return name;}
    
    /**
     * Accessor for symbol of unit
     */
    public String symbol() {return symbol;};

}