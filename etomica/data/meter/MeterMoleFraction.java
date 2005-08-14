package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Species;
import etomica.data.DataSourceScalar;
import etomica.data.Meter;
import etomica.units.Dimension;

/**
 * Meter for measurement of the species mole fraction in a phase.
 * Mole fraction is defined (number of molecules of species)/(number of molecules in phase).
 *
 * @author David Kofke
 */
public class MeterMoleFraction extends DataSourceScalar implements Meter {
    private Species species;
   
    public MeterMoleFraction() {
        super("Mole Fraction",Dimension.FRACTION);
    }
    
    public void setSpecies(Species s) {
        species = s;
    }
    public Species getSpecies() {return species;}

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species mole fraction in a phase");
        return info;
    }

    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
    	return (species == null) ? Double.NaN :
         	(double)phase.getAgent(species).moleculeCount()/(double)phase.moleculeCount();
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
}//end of MeterMoleFraction
