package etomica;

import etomica.units.Dimension;

/**
 * Meter for recording the total number of molecules in the phase
 */
public class MeterNMolecules extends MeterScalar implements EtomicaElement {
    
    private Species species;
    
    public MeterNMolecules() {
        this(Simulation.instance);
    }
    public MeterNMolecules(Simulation sim) {
        super(sim);
        setLabel("Molecules");
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Number of molecules in a phase");
        return info;
    }

    public void setSpecies(Species s) {species = s;}
    public Species getSpecies() {return species;}

    public Dimension getDimension() {return Dimension.QUANTITY;}

    public double getDataAsScalar(Phase phase) {
        return (species == null) ? phase.moleculeCount(): phase.getAgent(species).moleculeCount();
    }
 }
