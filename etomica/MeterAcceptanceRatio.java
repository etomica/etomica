package etomica;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.units.Count;

/**
 * Meter that reports the acceptance rate of MCMove type moves.
 * Returns acceptance rate as kept by the MCMove.
 */

public final class MeterAcceptanceRatio extends MeterAbstract {
    
    private MCMove move;
    
    public MeterAcceptanceRatio() {
        this(Simulation.instance);
    }
    public MeterAcceptanceRatio(Simulation sim) {
        super(sim, 1);
        setLabel("AcceptanceRatio");
    }
    public MeterAcceptanceRatio(Simulation sim, MCMove move) {
        this(sim);
        this.move = move;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Records the acceptance rate performed by the MCMove");
        return info;
    }
    public Unit defaultIOUnit() {return Count.UNIT;}
        
    /**
     * Returns dimensions of this meter's output, which in this case is QUANTITY.
     */
    public Dimension getDimension() {return Dimension.QUANTITY;}
        
    public void setMove(MCMove mv) {move = mv;}
        
    public MCMove getMove() {return move;}
        
       
    public void doMeasurement() {
        data[0] = move.acceptanceRatio();
    }    
}