package etomica;

import etomica.units.Dimension;

/**
 * Meter for recording the total number of molecules in the phase
 */
public class MeterNMolecules extends MeterAbstract implements EtomicaElement, MeterAtomic {
    
    private Species species;
    private SpeciesAgent speciesAgent;
    
    public MeterNMolecules() {
        this(Simulation.instance);
    }
    public MeterNMolecules(Simulation sim) {
        super(sim, 1);
        setLabel("Molecules");
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Number of molecules in a phase");
        return info;
    }

    public void setSpecies(Species s) {
        species = s;
        speciesAgent = species.getAgent(getPhase());
    }
    public Species getSpecies() {return species;}

    public Dimension getDimension() {return Dimension.QUANTITY;}

    public void doMeasurement() {
        if(speciesAgent == null) data[0] = phase.moleculeCount();
        else data[0] = speciesAgent.moleculeCount();
    }

    public double currentValue(Atom a) {
        return (a.node.parentSpecies()==species) ? 1.0 : 0.0;
    }

}
