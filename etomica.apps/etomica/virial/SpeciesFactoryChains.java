package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
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