package etomica.virial;

import etomica.chem.elements.ElementSimple;
import etomica.config.Conformation;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.space.Space;
import etomica.species.SpeciesSpheres;

/**
 * SpeciesFactory that makes a tangent sphere species.
 */
public class SpeciesFactoryTangentSpheres implements SpeciesFactory, java.io.Serializable {
    public SpeciesFactoryTangentSpheres(int nA, Conformation conformation) {
        this.nA = nA;
        this.conformation = conformation;
    }
    
    public ISpecies makeSpecies(ISimulation sim, Space space) {
        return new SpeciesSpheres(sim, nA, new ElementSimple(
                (sim.getSpeciesManager()).makeUniqueElementSymbol("TS"), 
                1.0), conformation, space);
    }
    
    private static final long serialVersionUID = 1L;
    private final int nA;
    private final Conformation conformation;
}
