/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;
import etomica.integrator.mcmove.MCMove;
import etomica.units.Decimal;
import etomica.units.dimensions.Fraction;
import etomica.units.Unit;

/**
 * Data source giving the acceptance rate of MCMove type moves.
 * Returns acceptance rate as kept by the MCMove.
 */

public class DataSourceAcceptanceRatio extends DataSourceScalar {
    
    protected MCMove move;

    /**
     * Constructs with no MCMove specified.  Calls to getDataAsScalar will return
     * Double.NaN until an MCMove is specified.
     */
    public DataSourceAcceptanceRatio() {
        this(null);
    }
    
    public DataSourceAcceptanceRatio(MCMove move) {
        super("AcceptanceRatio", Fraction.DIMENSION);
        setMove(move);
    }

    /**
     * @return Decimal.UNIT
     */
    public Unit defaultIOUnit() {return Decimal.UNIT;}
    
    /**
     * Sets the move for which the acceptance ratio is reported.
     */
    public void setMove(MCMove mv) {
        move = mv;
    }
        
    /**
     * @return the move for which the acceptance ratio is reported.
     */
    public MCMove getMove() {
        return move;
    }

    /**
     * Returns the value of the acceptance ratio as currently kept by
     * the MCMove class identified with setMove.
     */
    public double getDataAsScalar() {
        if (move == null) return Double.NaN;
        return move.getTracker().acceptanceRatio();
    }
    
    private static final long serialVersionUID = 1L;
}
