package etomica.virial;

import etomica.chem.elements.ElementSimple;
import etomica.config.Conformation;
import etomica.simulation.ISimulation;
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
    
    public Species makeSpecies(ISimulation sim) {
        return new SpeciesSpheres(sim, nA, new ElementSimple(
                (sim.getSpeciesManager()).makeUniqueElementSymbol("TS"), 
                1.0), conformation);
    }
    
    private static final long serialVersionUID = 1L;
    private final int nA;
    private final Conformation conformation;
}
