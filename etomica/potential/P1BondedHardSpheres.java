package etomica.potential;

import etomica.simulation.Simulation;

/**
 * Intramolecular potential in which bonded and nonbonded atoms interact with a
 * hard potential (P2HardBond and P2HardSphere by default).
 * 
 * @author David Kofke
 */
 
public class P1BondedHardSpheres {

    public static PotentialGroup makeP1BondedHardSpheres(Simulation sim) {
        return P1IntraSimple.makeP1IntraSimple(sim.getPotentialMaster(), new P2HardBond(sim), new P2HardSphere(sim));
    }
}
