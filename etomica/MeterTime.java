package etomica;
import etomica.units.Dimension;

/**
 * Meter for the simulation time, as given by an IntegratorMD instance.
 * 
 * @author David Kofke
 */
 
 /* History of changes
  * 7/29/02 (DAK) new
  */

public class MeterTime extends MeterAbstract implements EtomicaElement {
        
    double time0;
    
    public MeterTime() {
        this(Simulation.instance);
    }
    public MeterTime(Simulation sim) {
        super(sim, 1);
        this.reset();
        setLabel("Elapsed time");
    }
                
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Reports elapsed time in a molecular dynamics simulation");
        return info;
    }

    /**
     * Indicates units (time) of measured property.
     */
    public Dimension getDimension() {return Dimension.TIME;}
        
    /**
     * Sets the meter's current time as the new zero for subsequent values.
     */
    public void reset() {time0 = (integrator == null) ? 0.0 : ((IntegratorMD)integrator).elapsedTime();}
        
    /**
     * Overrides parent class method to cause no action to be performed in response 
     * to interval event.
     */
    public void intervalAction(Integrator.IntervalEvent evt) {}
    /**
     * Returns the simulation time elapsed since the instantiation of
     * this meter, or since the last call to reset().
     */
    public void doMeasurement() {
    	data[0] = ((IntegratorMD)integrator).elapsedTime() - time0;
    }

}//end of MeterTime
