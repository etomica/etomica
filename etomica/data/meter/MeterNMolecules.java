package etomica.data.meter;

import etomica.Atom;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Species;
import etomica.data.Data;
import etomica.data.DataSourceAtomic;
import etomica.data.DataSourceScalar;
import etomica.data.types.DataDouble;
import etomica.units.Dimension;

/**
 * Meter for recording the total number of molecules in the phase
 */
public class MeterNMolecules extends DataSourceScalar implements DataSourceAtomic, Meter {
    
    private Species species;
    
    public MeterNMolecules() {
        super("Molecules",Dimension.QUANTITY);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Number of molecules in a phase");
        return info;
    }

    public Class getDataType() {
        return DataDouble.class;
    }
    
    public void setSpecies(Species s) {species = s;}
    public Species getSpecies() {return species;}

    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
        return (species == null) ? phase.moleculeCount(): phase.getAgent(species).moleculeCount();
    }
    
    public Data getData(Atom atom) {
        data.x = (species == null || atom.type.getSpecies() == species) ? 1 : 0;
        return data;
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
}
