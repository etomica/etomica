package etomica.virial;

import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;


/**
 * SpeciesFactory that makes SpeciesSpheres 
 */
public class SpeciesFactoryChains implements SpeciesFactory {
    
    public Species makeSpecies(Simulation sim) {
        return new SpeciesSpheres(sim);
    }
    
    public int nA;
}