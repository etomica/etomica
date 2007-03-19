package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.data.DataSourceScalar;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.species.Species;
import etomica.units.DimensionRatio;
import etomica.units.Quantity;
import etomica.units.Volume;

/**
 * Meter for measurement of the total molecule number density in a phase
 * Molecule number density is defined (number of molecules)/(volume of phase)
 */
public class MeterDensity extends DataSourceScalar {
    
    public MeterDensity(Space space) {
        super("Number Density",new DimensionRatio(Quantity.DIMENSION, Volume.dimension(space.D())));
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
        			phase.getAgent(species).getNMolecules())
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

    private static final long serialVersionUID = 1L;
    private Phase phase;
    private Species species;
}
