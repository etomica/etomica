package etomica;

/**
 * Meter for measurement of the species mole fraction in a phase.
 * Mole fraction is defined (number of molecules of species)/(number of molecules in phase).
 *
 * @author David Kofke
 */
public class MeterMoleFraction extends MeterScalar implements EtomicaElement {
    private Species species;
   
    public MeterMoleFraction() {
        super();
        setLabel("Mole fraction");
    }
    
    public void setSpecies(Species s) {
        species = s;
    }
    public Species getSpecies() {return species;}

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species mole fraction in a phase");
        return info;
    }

    public double getDataAsScalar(Phase phase) {
    	return (species == null) ? Double.NaN :
         	(double)phase.getAgent(species).moleculeCount()/(double)phase.moleculeCount();
     }
        
    public etomica.units.Dimension getDimension() {return etomica.units.Dimension.FRACTION;}

}//end of MeterMoleFraction
