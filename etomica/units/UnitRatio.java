package etomica.units;

/**
 * Class for constructing units forms as ratios of other units. 
 * Examples include "Joules per mole", "gram per liter", and so on.
 * Conversion factor and labels are determined from the values given by the
 * constituent units.
 */
public class UnitRatio extends BaseUnit {
        
    /**
     * Constructs the compound unit using conversion factors and labels from the
     * constituent units.  These factors/labels are computed once during construction,
     * and are not updated if the units are modified by changing their prefix value.
     * @param numerator the unit in the numerator of the compound unit
     * @param denominator the unit in the denominator of the compount unit
     */
    public UnitRatio(Unit numerator, Unit denominator) {
        super(numerator.toSim(1.0)/denominator.toSim(1.0),                        //baseTo
              numerator.toString() + " per " + denominator.toString(),            //name
              numerator.symbol() + "/" + denominator.symbol(),                    //symbol
			  new DimensionRatio(numerator.dimension(), denominator.dimension()), //dimension
              Prefix.ALLOWED                                                 //prefixAllowed
              );
    }
    /**
     * @param numerator the unit in the numerator of the compound unit
     * @param denominator the unit in the denominator of the compount unit
     * @param name a description of the unit (e.g., Joules per mole)
     * @param symbol an abbreviated description of the unit (e.g., J/mol)
     */
    public UnitRatio(Unit numerator, Unit denominator, String name, String symbol) {
		super(numerator.toSim(1.0)/denominator.toSim(1.0),                        //baseTo
			  name,            													  //name
			  symbol,                    										  //symbol
			  new DimensionRatio(numerator.dimension(), denominator.dimension()), //dimension
			  Prefix.ALLOWED                                                 //prefixAllowed
			  );
    }
}
    
