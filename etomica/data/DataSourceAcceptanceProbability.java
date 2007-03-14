package etomica.data;

import etomica.integrator.mcmove.MCMove;
import etomica.units.Decimal;
import etomica.units.Unit;

/**
 * Data source giving the average acceptance probability of MCMove type moves.
 * Returns acceptance probability as kept by the MCMove.
 */

public class DataSourceAcceptanceProbability extends DataSourceAcceptanceRatio {

    /**
     * Constructs with no MCMove specified.  Calls to getDataAsScalar will return
     * Double.NaN until an MCMove is specified.
     *
     */
    public DataSourceAcceptanceProbability() {
        this(null);
    }
    public DataSourceAcceptanceProbability(MCMove move) {
        super(move);
    }
    
    /**
     * @return Decimal.UNIT
     */
    public Unit defaultIOUnit() {return Decimal.UNIT;}
    
   /**
     * Returns the value of the acceptance ratios as currently kept by
     * the MCMove classes identified with setMove.  Each element of the
     * returned array applies to the corresponding element in the
     * MCMove array.
     */
    public double getDataAsScalar() {
        if (move == null) return Double.NaN;
        return move.getTracker().acceptanceProbability();
    }

    private static final long serialVersionUID = 1L;
}
