package etomica.potential;

import etomica.EtomicaInfo;
import etomica.simulation.Simulation;

/**
 * Intramolecular potential in which bonded atoms interact with a hard tether
 * potential, and nonbonded atoms interact as hard spheres.
 *
 * @author David Kofke
 */
 
public class P1TetheredHardSpheres extends P1IntraSimple {
    
    public P1TetheredHardSpheres(Simulation sim) {
        super(sim.space, sim.potentialMaster, new P2HardSphere(sim), new P2Tether(sim));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Bonded atoms tethered, nonbonded interact as hard spheres");
        return info;
    }
}
   
