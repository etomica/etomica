package etomica.data;
import etomica.EtomicaInfo;
import etomica.integrator.mcmove.MCMove;
import etomica.units.Decimal;
import etomica.units.Fraction;
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
   
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Records the acceptance rate performed by the MCMove");
        return info;
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