package etomica.potential;

import etomica.simulation.Simulation;

/**
 * Intramolecular potential in which bonded atoms interact with a hard tether
 * potential, and nonbonded atoms interact as hard spheres.
 *
 * @author David Kofke
 */
 
public class P1TetheredHardSpheres {
    
    public PotentialGroup makeP1TetheredHardSpheres(Simulation sim) {
        return P1IntraSimple.makeP1IntraSimple(sim.getPotentialMaster(), new P2HardSphere(sim), new P2Tether(sim));
    }
}
   
