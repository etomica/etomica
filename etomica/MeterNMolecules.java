package etomica;

import etomica.units.Dimension;

/**
 * Meter for recording the total number of molecules in the phase
 */
public class MeterNMolecules extends Meter
{
    public MeterNMolecules() {
        this(Simulation.instance);
    }
    public MeterNMolecules(Simulation sim)
    {
        super(sim);
        setLabel("Molecules");
    }
    
    /**
     * Declaration that this meter does not use the boundary object of phase when making its measurements
     */
    public final boolean usesPhaseBoundary() {return false;}
    /**
     * Declaration that this meter does not use the iteratorFactory of phase when making its measurements
     */
    public final boolean usesPhaseIteratorFactory() {return false;}

    public Dimension getDimension() {return Dimension.QUANTITY;}

    public double currentValue()
    {
        return phase.moleculeCount;
    }

}
