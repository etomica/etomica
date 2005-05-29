package etomica.data;

import etomica.integrator.MCMove;

/**
 * Data source giving the average acceptance probability of MCMove type moves.
 * Returns acceptance probability as kept by the MCMove.
 */

public class DataSourceAcceptanceProbability extends DataSourceAcceptanceRatio {

    public DataSourceAcceptanceProbability() {
        this(new MCMove[0]);
    }
    public DataSourceAcceptanceProbability(MCMove[] move) {
        super(move);
        setLabel("AcceptanceProbability");
    }
    
    /**
     * Returns the value of the acceptance ratios as currently kept by
     * the MCMove classes identified with setMove.  Each element of the
     * returned array applies to the corresponding element in the
     * MCMove array.
     */
    public double[] getData() {
        for (int i=0; i<move.length; i++) {
            ratioArray[i] = move[i].acceptanceProbability();
        }
        return ratioArray;
    }
}
