package etomica;//
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.units.Count;

/**
 * Meter that keeps track of the acceptance rate of MCMove type moves
 * 
 * Methods average and error are meaningless for this integrator, and return Not-a-Number
 */

public final class MeterAcceptanceRatio extends MeterScalar {
    
    private MCMove move;
    
    public MeterAcceptanceRatio() {
        this(Simulation.instance);
    }
    public MeterAcceptanceRatio(Simulation sim) {
        super(sim);
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
    public Unit defaultIOUnit() {return new Unit(Count.UNIT);}
        
    /**
     * Returns dimensions of this meter's output, which in this case is QUANTITY.
     */
    public Dimension getDimension() {return Dimension.QUANTITY;}
    
    
    /**
     * Returns Not-a-Number
     */
    public double average() {return Double.NaN;}
    /**
     * Returns Not-a-Number
     */
    public double error() {return Double.NaN;}
    
    /**
     * Returns currentValue
    */ 
    public double mostRecent() {return currentValue();}
    
    
    public void setMove(MCMove mv) {move = mv;}
        
    public MCMove getMove() {return move;}
        
       
    public double currentValue() {
        return move.acceptanceRatio();
    }    
}