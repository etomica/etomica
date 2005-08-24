/*
 * Created on May 18, 2005
 */
package etomica.models.hexane;

import etomica.atom.AtomSequencerFactory;
import etomica.simulation.Simulation;

/**
 * Species used to create hexane molecules per Dr. Monson's data. Hydrogen
 * molecules are ignored.
 * 
 * @author nancycribbin
 */

public class SpeciesHexane extends etomica.species.SpeciesSpheres {
    public SpeciesHexane(Simulation sim){
        super(sim, sim.potentialMaster.sequencerFactory(), 6, 
                new ConformationHexane(sim.space));
    }
    
    public SpeciesHexane(Simulation sim, AtomSequencerFactory seqFactory) {
        super(sim, seqFactory, 6, new ConformationHexane(sim.space));
    }
    
}
