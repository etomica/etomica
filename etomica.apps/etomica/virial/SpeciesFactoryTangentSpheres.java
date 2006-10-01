package etomica.virial;

import etomica.atom.AtomTypeRoot;
import etomica.chem.elements.ElementSimple;
import etomica.config.Conformation;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryTangentSpheres implements SpeciesFactory, java.io.Serializable {
    public SpeciesFactoryTangentSpheres(int nA, Conformation conformation) {
        this.nA = nA;
        this.conformation = conformation;
    }
    
    public Species makeSpecies(Simulation sim) {
        return new SpeciesSpheres(sim, nA, new ElementSimple(
                ((AtomTypeRoot)sim.speciesRoot.type).makeUniqueElementSymbol("TS"), 
                sim.getDefaults().atomMass), conformation);
    }
    
    private final int nA;
    private final Conformation conformation;
}
