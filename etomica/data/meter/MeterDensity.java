package etomica.data.meter;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Species;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;

/**
 * Meter for measurement of the total molecule number density in a phase
 * Molecule number density is defined (number of molecules)/(volume of phase)
 *
 * MIGHT WANT TO CHANGE THIS TO A METER.RATIO
 */
public class MeterDensity extends MeterScalar implements EtomicaElement
{
    public MeterDensity() {
        super();
        setLabel("Number density");
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Number density (molecules/volume) in a phase");
        return info;
    }
    
    public void setSpecies(Species s) {
        species = s;
    }
    public Species getSpecies() {
    	return species;
    }

    public double getDataAsScalar(Phase phase) {
        return (species == null ? 
        			phase.moleculeCount() : 
        			phase.getAgent(species).moleculeCount())
				/phase.volume();
    }
    
    public Dimension getDimension() {return new DimensionRatio(Dimension.QUANTITY, Dimension.VOLUME);}

    private Species species;
}
