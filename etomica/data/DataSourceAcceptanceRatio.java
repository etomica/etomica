package etomica.data;
import etomica.EtomicaInfo;
import etomica.integrator.MCMove;
import etomica.units.Decimal;
import etomica.units.Dimension;
import etomica.units.Unit;

/**
 * Data source giving the acceptance rate of MCMove type moves.
 * Returns acceptance rate as kept by the MCMove.
 */

public class DataSourceAcceptanceRatio extends DataSourceAdapter {
    
    protected MCMove[] move;
    protected double[] ratioArray;
    
    public DataSourceAcceptanceRatio() {
    	this(new MCMove[0]);
    }
    public DataSourceAcceptanceRatio(MCMove[] move) {
        super(Dimension.FRACTION);
        setLabel("AcceptanceRatio");
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
     * Sets the moves for which the acceptance ratios are reported.
     * @param mv
     */
    public void setMove(MCMove[] mv) {
        move = (MCMove[])mv.clone();
        ratioArray = new double[move.length];
    }
        
    /**
     * @return the moves for which the acceptance ratios are reported.
     */
    public MCMove[] getMove() {return move;}

    /**
     * Returns the value of the acceptance ratios as currently kept by
     * the MCMove classes identified with setMove.  Each element of the
     * returned array applies to the corresponding element in the
     * MCMove array.
     */
    public double[] getData() {
    	for (int i=0; i<move.length; i++) {
    		ratioArray[i] = move[i].acceptanceRatio();
    	}
        return ratioArray;
    }
    
    public int getDataLength() {
        return move.length;
    }
}