package etomica.data.meter;

import etomica.DataInfo;
import etomica.EtomicaInfo;
import etomica.Meter;
import etomica.Phase;
import etomica.Species;
import etomica.data.DataSourceScalar;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;

/**
 * Meter for measurement of the total molecule number density in a phase
 * Molecule number density is defined (number of molecules)/(volume of phase)
 */
public class MeterDensity extends DataSourceScalar implements Meter {
    
    public MeterDensity() {
        super(new DataInfo("Number Density",new DimensionRatio(Dimension.QUANTITY,Dimension.VOLUME)));
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

    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
        return (species == null ? 
        			phase.moleculeCount() : 
        			phase.getAgent(species).moleculeCount())
				/phase.volume();
    }
    
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
    private Species species;
}
