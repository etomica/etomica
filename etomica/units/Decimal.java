/*
 * History
 * Created on Jul 30, 2004 by kofke
 */
package etomica.units;

/**
 * @author kofke
 *
 * Decimal representation of something that represents the fractional 
 * amount of a whole (e.g., mole fraction) as a decimal value typically
 * between 0 and 1.
 */
public class Decimal extends etomica.units.BaseUnit.Fraction {

  /**
   * Singleton instance of this unit 
   */
	public static final Decimal UNIT = new Decimal();

	private Decimal() {
       super(
        	1.0,
        	"Decimal",
        	"",
        	Prefix.NOT_ALLOWED
        	);
	}
}
