package etomica.units;

/**
 * Class for constructing units forms as ratios of other units. 
 * Examples include "Joules per mole", "gram per liter", and so on.
 * Conversion factor and labels are determined from the values given by the
 * constituent units.
 */
public class UnitRatio extends Unit {
        
    private double toPixels;

    /**
     * Constructs the compound unit using conversion factors and labels from the
     * constituent units.  These factors/labels are computed once during construction,
     * and are not updated if the units are modified by changing their prefix value.
     * @param numerator the unit in the numerator of the compound unit
     * @param denominator the unit in the denominator of the compount unit
     */
    public UnitRatio(Unit numerator, Unit denominator) {
        super(Prefix.NULL,                                                        //prefix
              numerator.toSim(1.0)/denominator.toSim(1.0),                        //baseTo
              numerator.toString() + " per " + denominator.toString(),            //name
              numerator.symbol() + "/" + denominator.symbol(),                    //symbol
              Prefix.NOT_ALLOWED,                                                 //prefixAllowed
              new DimensionRatio(numerator.dimension(), denominator.dimension()) //dimension
              );
    }
    /**
     * @param numerator the unit in the numerator of the compound unit
     * @param denominator the unit in the denominator of the compount unit
     * @param name a description of the unit (e.g., Joules per mole)
     * @param symbol an abbreviated description of the unit (e.g., J/mol)s
     */
    public UnitRatio(Unit numerator, Unit denominator, String name, String symbol) {
        this(numerator, denominator);
        this.name = name;
        this.symbol = symbol;
    }
    
    public void setToPixels(double t) {toPixels = t;}
    public double toPixels(double x) {return x*toPixels;}
}