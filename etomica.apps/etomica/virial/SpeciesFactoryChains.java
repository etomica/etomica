package etomica.virial;

import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.space.Space;
import etomica.species.SpeciesSpheres;


/**
 * SpeciesFactory that makes SpeciesSpheres 
 */
public class SpeciesFactoryChains implements SpeciesFactory {
    
    public ISpecies makeSpecies(ISimulation sim, Space _space) {
        return new SpeciesSpheres(sim, _space);
    }
    
    public int nA;
}