package etomica;

import etomica.units.*;

/**
 * Meter for measurement of the total molecule number density in a phase
 * Molecule number density is defined (number of molecules)/(volume of phase)
 *
 * MIGHT WANT TO CHANGE THIS TO A METER.RATIO
 */
public class MeterDensity extends MeterScalar implements EtomicaElement
{
    public MeterDensity() {
        this(Simulation.instance);
    }
    public MeterDensity(Simulation sim) {
        super(sim);
        setLabel("Number density");
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Number density (molecules/volume) in a phase");
        return info;
    }

    public double getDataAsScalar(Phase phase) {
        return phase.moleculeCount()/phase.volume();
   }
    
    public Dimension getDimension() {return new DimensionRatio(Dimension.QUANTITY, Dimension.VOLUME);}
}
