package etomica.virial;

import etomica.simulation.ISimulation;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheres;


/**
 * SpeciesFactory that makes SpeciesSpheres 
 */
public class SpeciesFactoryChains implements SpeciesFactory {
    
    public ISpecies makeSpecies(ISimulation sim) {
        return new SpeciesSpheres(sim);
    }
    
    public int nA;
}