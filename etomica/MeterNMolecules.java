package etomica;

import etomica.units.Dimension;

/**
 * Meter for recording the total number of molecules in the phase
 */
public class MeterNMolecules extends Meter implements EtomicaElement, Meter.Atomic {
    
    private Species species;
    private Species.Agent speciesAgent;
    
    public MeterNMolecules() {
        this(Simulation.instance);
    }
    public MeterNMolecules(Simulation sim)
    {
        super(sim);
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

    /**
     * Declaration that this meter does not use the boundary object of phase when making its measurements
     */
    public final boolean usesPhaseBoundary() {return false;}
    /**
     * Declaration that this meter does not use the iteratorFactory of phase when making its measurements
     */
    public final boolean usesPhaseIteratorFactory() {return false;}

    public Dimension getDimension() {return Dimension.QUANTITY;}

    public double currentValue() {
        if(speciesAgent == null) return phase.moleculeCount;
        else return speciesAgent.moleculeCount();
    }

    public double currentValue(Atom a) {
        return (a.parentMolecule().parentSpecies()==species) ? 1.0 : 0.0;
    }

}
