package etomica;

import etomica.units.*;

/**
 * Meter for measurement of the total molecule number density in a phase
 * Molecule number density is defined (number of molecules)/(volume of phase)
 *
 * MIGHT WANT TO CHANGE THIS TO A METER.RATIO
 */
public class MeterDensity extends Meter implements Meter.Atomic, EtomicaElement
{
    public MeterDensity() {
        this(Simulation.instance);
    }
    public MeterDensity(Simulation sim)
    {
        super(sim);
        setLabel("Number density");
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Number density (molecules/volume) in a phase");
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
        return phase.moleculeCount()/phase.volume();
    }
    
    public double currentValue(Atom a) {
        return 1/phase.volume();
    }

    public Dimension getDimension() {return new DimensionRatio(Dimension.QUANTITY, Dimension.VOLUME);}
}
