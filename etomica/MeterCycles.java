package etomica;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.units.Count;

/**
 * Meter that keeps track of the number of intervalEvents fired by an integrator.
 * More precisely, currentValue method returns number of events received divided by updateInterval, 
 * which is set to 1 by default.
 * Methods average and error are meaningless for this integrator, and return Not-a-Number
 */

public final class MeterCycles extends MeterScalar implements EtomicaElement
{
    private int count = 0;
    
    public MeterCycles() {
        this(Simulation.instance);
    }
    public MeterCycles(Simulation sim){
        super(sim);
        setLabel("Cycles");
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Records the number of simulation cycles performed by the integrator");
        return info;
    }

    public Unit defaultIOUnit() {return Count.UNIT;}
        
    /**
     * Returns dimensions of this meter's output, which in this case is QUANTITY.
     */
    public Dimension getDimension() {return Dimension.QUANTITY;}
    
    /**
     * Resets counter to zero
     */
    public void reset() {count = 0;}
    /**
     * Returns the number of interval events received, divided by updateInterval
     */
    public double getDataAsScalar() {
    	count += getUpdateInterval();
        return (double)count;  
    }
}