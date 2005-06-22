package etomica.data;

import etomica.integrator.MCMove;

/**
 * Data source giving the average acceptance probability of MCMove type moves.
 * Returns acceptance probability as kept by the MCMove.
 */

public class DataSourceAcceptanceProbability extends DataSourceAcceptanceRatio {

    public DataSourceAcceptanceProbability() {
        this(null);
    }
    public DataSourceAcceptanceProbability(MCMove move) {
        super(move);
    }
    
    /**
     * Returns the value of the acceptance ratios as currently kept by
     * the MCMove classes identified with setMove.  Each element of the
     * returned array applies to the corresponding element in the
     * MCMove array.
     */
    public double getDataAsScalar() {
        if (move == null) return Double.NaN;
        return move.acceptanceProbability();
    }
}
