package etomica;
import etomica.units.Dimension;

/**
 * Meter for measurement of the temperature based on kinetic-energy equipartition
 */

public final class MeterTemperature extends Meter
{
    public MeterTemperature() {
        this(Simulation.instance);
    }
    public MeterTemperature(Simulation sim)
    {
        super(sim);
        setLabel("Temperature");
    }

    /**
     * Declaration that this meter does not use the boundary object of phase when making its measurements
     */
    public final boolean usesPhaseBoundary() {return false;}
    /**
     * Declaration that this meter does not use the iteratorFactory of phase when making its measurements
     */
    public final boolean usesPhaseIteratorFactory() {return false;}
        
    public double currentValue()
    {
        return value(phase);
    }
    
	public Dimension getDimension() {return Dimension.TEMPERATURE;}

/**
 * Class method to compute the temperature of a phase from its total kinetic energy using equipartition
 */
    public static double value(Phase phase) {
        return (2./(double)(phase.atomCount*phase.parentSimulation().space().D()))*phase.energy.kinetic();
    }    
}