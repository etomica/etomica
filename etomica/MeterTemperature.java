package etomica;
import etomica.units.Dimension;

/**
 * Meter for measurement of the temperature based on kinetic-energy equipartition
 */

public final class MeterTemperature extends Meter implements EtomicaElement
{
    public MeterTemperature() {
        this(Simulation.instance);
    }
    public MeterTemperature(Simulation sim)
    {
        super(sim);
        setLabel("Temperature");
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Records temperature as given via kinetic energy");
        return info;
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
        return currentValue(phase);
    }
    
	public Dimension getDimension() {return Dimension.TEMPERATURE;}

/**
 * Class method to compute the temperature of a phase from its total kinetic energy using equipartition
 */
    public static double currentValue(Phase phase) {
        return (2./(double)(phase.atomCount()*phase.parentSimulation().space().D()))*MeterKineticEnergy.currentValue(phase);
    }    
}