/*
 * History
 * Created on Jul 30, 2004 by kofke
 */
package etomica.units;

/**
 * @author kofke
 *
 * Decimal representation of something that represents the fractional 
 * amount of a whole (e.g., mole fraction) as a percentage value typically
 * between 0 and 100.
 */
public class Percent extends etomica.units.BaseUnit.Fraction {

  /**
   * Singleton instance of this unit 
   */
	public static final Percent UNIT = new Percent();

	private Percent() {
       super(
        	0.01,
        	"Percent",
        	"%",
        	Prefix.NOT_ALLOWED
        	);
	}
}
