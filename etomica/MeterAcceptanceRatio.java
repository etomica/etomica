package etomica;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.units.Count;

/**
 * Meter that reports the acceptance rate of MCMove type moves.
 * Returns acceptance rate as kept by the MCMove.
 */

public final class MeterAcceptanceRatio extends DataSourceAdapter {
    
    private MCMove[] move = new MCMove[0];
    private String label;
    private double[] ratioArray;
    
    public MeterAcceptanceRatio() {
    	super(Dimension.FRACTION);
        setLabel("AcceptanceRatio");
    }
    public MeterAcceptanceRatio(Simulation sim, MCMove[] move) {
        this();
        setMove(move);
    }

    public void setLabel(String label) {
    	this.label = label;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Records the acceptance rate performed by the MCMove");
        return info;
    }
    public Unit defaultIOUnit() {return Count.UNIT;}
        
    public void setMove(MCMove[] mv) {
        move = (MCMove[])mv.clone();
        ratioArray = new double[move.length];
    }
        
    public MCMove[] getMove() {return move;}

    public double[] getData() {
    	for (int i=0; i<move.length; i++) {
    		ratioArray[i] = move[i].acceptanceRatio();
    	}
        return ratioArray;
    }
}