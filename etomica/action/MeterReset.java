package etomica.action;

import etomica.*;

/**
 * Calls the reset method on the given meters, which rezeros 
 * all accumulators used to report averages and confidence limits.
 */
public final class MeterReset extends MeterAction {
    
    public MeterReset() {
        this(Simulation.instance);
    }
    /**
     * Defines target meters via the given array.
     */
    public MeterReset(MeterAbstract[] meters) {
        super(meters);
        setLabel("Reset averages");
    }
    /**
     * Defines target meters as those in the meterList of the given simulation.
     */
    public MeterReset(Simulation sim) {
        super(sim);
        setLabel("Reset averages");
    }
    
    /**
     * Resets the meters in the given array.
     */
    public static void doAction(MeterAbstract[] meters) {
        for(int i=0; i<meters.length; i++) {
            doAction(meters[i]);
        }
    }
    /**
     * Resets the given meter.
     */
    public static void doAction(MeterAbstract meter) {
        meter.reset();
    }
    
    /**
     * Resets the given meter.  Implementation of abstract Action method.
     * Normally called as an action listener.
     */
    public void actionPerformed(MeterAbstract meter) {
        MeterReset.doAction(meter);
    }
}
        