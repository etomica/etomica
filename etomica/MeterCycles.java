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

public final class MeterCycles extends Meter implements EtomicaElement
{
    private int count = 0;
    
    public MeterCycles() {
        this(Simulation.instance);
    }
    public MeterCycles(Simulation sim)
    {
        super(sim);
        setLabel("Cycles");
        setUpdateInterval(1);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Records the number of simulation cycles performed by the integrator");
        return info;
    }

    public Unit defaultIOUnit() {return new Unit(Count.UNIT);}
        
    /**
     * Returns dimensions of this meter's output, which in this case is QUANTITY.
     */
    public Dimension getDimension() {return Dimension.QUANTITY;}
    
    /**
     * Increments the counter.  Normally called by intervalAction method defined in superclass.
     */
    public void updateSums() {count += getUpdateInterval();}
    /**
     * Resets counter to zero
     */
    public void reset() {count = 0;}
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
    
    /**
     * Returns the number of interval events received, divided by updateInterval
     */
    public double currentValue()
    {
        return (double)count;  
    }
}