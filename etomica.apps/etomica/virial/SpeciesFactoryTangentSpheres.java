package etomica.virial;

import etomica.chem.elements.ElementSimple;
import etomica.api.IConformation;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.simulation.SpeciesManager;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheres;

/**
 * SpeciesFactory that makes a tangent sphere species.
 */
public class SpeciesFactoryTangentSpheres implements SpeciesFactory, java.io.Serializable {
    public SpeciesFactoryTangentSpheres(int nA, IConformation conformation) {
        this.nA = nA;
        this.conformation = conformation;
    }
    
    public ISpecies makeSpecies(ISimulation sim, ISpace space) {
        return new SpeciesSpheres(sim, nA, new ElementSimple(
                ((SpeciesManager)sim.getSpeciesManager()).makeUniqueElementSymbol("TS"), 
                1.0), conformation, space);
    }
    
    private static final long serialVersionUID = 1L;
    private final int nA;
    private final IConformation conformation;
}
