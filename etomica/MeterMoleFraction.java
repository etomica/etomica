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
        this(Simulation.instance);
    }
    public MeterMoleFraction(Simulation sim) {
        super(sim);
        setLabel("Mole fraction");
    }
    
    public void setSpecies(Species s) {
        species = s;
//        speciesAgent = species.getAgent(getPhase());
    }
    public Species getSpecies() {return species;}

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species mole fraction in a phase");
        return info;
    }

    public void getData(Phase phase) {
        if(species == null) phaseData[0] = Double.NaN;
        else {
        	SpeciesAgent agent = phase.getAgent(species);
        	phaseData[0] = (double)agent.moleculeCount()/(double)phase.moleculeCount();
        }
     }
    
    public double currentValue(Atom a) {
        return (a.node.parentSpecies()==species) ? 1.0/(double)phase.moleculeCount() : 0.0;
    }
    
    public etomica.units.Dimension getDimension() {return etomica.units.Dimension.NULL;}

}//end of MeterMoleFraction
