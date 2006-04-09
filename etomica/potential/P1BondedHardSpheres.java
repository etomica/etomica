package etomica.potential;

import etomica.EtomicaInfo;
import etomica.simulation.Simulation;

/**
 * Intramolecular potential in which bonded and nonbonded atoms interact with a
 * hard potential (P2HardBond and P2HardSphere by default).
 * 
 * @author David Kofke
 */
 
public class P1BondedHardSpheres extends P1IntraSimple {

    public P1BondedHardSpheres(Simulation sim) {
        super(sim.space, sim.potentialMaster, new P2HardBond(sim), new P2HardSphere(sim));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Bonded atoms tethered, nonbonded interact as hard spheres");
        info.getFeatures().add("NBODIES", 1);
        return info;
    }

}
