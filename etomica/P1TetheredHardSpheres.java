package etomica;

/**
 * Intramolecular potential in which bonded atoms interact with a hard tether
 * potential, and nonbonded atoms interact as hard spheres.
 *
 * @author David Kofke
 */
 
public class P1TetheredHardSpheres extends Potential1Group implements EtomicaElement, Potential1.Intramolecular {
    
    public String getVersion() {return "P1TetheredHardSpheres:01.11.05/"+Potential1Group.VERSION;}
    
    public final P2HardSphere p2HardSphere;
    public final P2Tether p2Tether;
    
    public P1TetheredHardSpheres() {
        this(Simulation.instance.hamiltonian.potential);
    }
    
    public P1TetheredHardSpheres(PotentialGroup parent) {
        super(parent);
        p2HardSphere = new P2HardSphere(this);
        p2Tether = new P2Tether(this);
	    p2Tether.setIterator(new ApiGeneral(parentSimulation().space,
	            new AtomIteratorSequential(false),
	            new AtomIteratorBonds()));
	    p2HardSphere.setIterator(new ApiGeneral(parentSimulation().space,
	            new AtomIteratorSequential(false),
	            new AtomIteratorNonbonded()));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Bonded atoms tethered, nonbonded interact as hard spheres");
        return info;
    }


    public double energy(Atom a) {
        return 0.0;
    }

    
}//end of P1TetheredHardSpheres
   
