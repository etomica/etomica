package etomica;

/**
 * Meter for measurement of the species molecule number density in a phase
 * Molecule number density is defined (number of molecules of species)/(volume of phase)
 *
 */
public class MeterDensitySpecies extends MeterDensity implements EtomicaElement
{
    private Species species;
    private SpeciesAgent speciesAgent;
    
    public MeterDensitySpecies() {
        this(Simulation.instance);
    }
    public MeterDensitySpecies(Simulation sim) {
        super(sim);
        setLabel("Number density");
    }
    public MeterDensitySpecies(Species s) {
        this(s.parentSimulation());
        setSpecies(s);
    }
    
    public void setSpecies(Species s) {
        species = s;
        speciesAgent = species.getAgent(getPhase());
    }
    public Species getSpecies() {return species;}

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Species number density (molecules/volume) in a phase");
        return info;
    }

    public double currentValue() {
        if(species == null) return 0.0;
        return speciesAgent.moleculeCount()/phase.volume();
    }
    
    public double currentValue(Atom a) {
        return (a.parentSpecies()==species) ? 1/phase.volume() : 0.0;
    }

}
