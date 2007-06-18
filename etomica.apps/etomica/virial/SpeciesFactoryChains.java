package etomica.virial;

import etomica.simulation.ISimulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;


/**
 * SpeciesFactory that makes SpeciesSpheres 
 */
public class SpeciesFactoryChains implements SpeciesFactory {
    
    public Species makeSpecies(ISimulation sim) {
        return new SpeciesSpheres(sim);
    }
    
    public int nA;
}