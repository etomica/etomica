package etomica.virial;

import etomica.chem.elements.ElementSimple;
import etomica.config.Conformation;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;

/**
 * SpeciesFactory that makes a tangent sphere species.
 */
public class SpeciesFactoryTangentSpheres implements SpeciesFactory, java.io.Serializable {
    public SpeciesFactoryTangentSpheres(int nA, Conformation conformation) {
        this.nA = nA;
        this.conformation = conformation;
    }
    
    public Species makeSpecies(Simulation sim) {
        return new SpeciesSpheres(sim, nA, new ElementSimple(
                (sim.getSpeciesManager()).makeUniqueElementSymbol("TS"), 
                sim.getDefaults().atomMass), conformation);
    }
    
    private static final long serialVersionUID = 1L;
    private final int nA;
    private final Conformation conformation;
}
