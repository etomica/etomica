package etomica;

/**
 * Intramolecular potential in which bonded atoms interact with a hard tether
 * potential, and nonbonded atoms interact as hard spheres.
 *
 * @author David Kofke
 */
 
public class P1TetheredHardSpheres extends PotentialGroup implements Potential1.Intramolecular {
    
    public String getVersion() {return "P1TetheredHardSpheres:01.11.05/"+PotentialGroup.VERSION;}
    
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
	            new AtomIteratorList(),
	            new AtomIteratorBonds()));
	    p2HardSphere.setIterator(new ApiGeneral(parentSimulation().space,
	            new AtomIteratorList(),
	            new AtomIteratorNonbonded(parent.parentSimulation())));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Bonded atoms tethered, nonbonded interact as hard spheres");
        return info;
    }

    public double energy(Atom a) {
        return 0.0;
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic();
        Simulation.instance = sim;
	    IntegratorHard integratorHard = new IntegratorHard();
	    SpeciesSpheres speciesSpheres = new SpeciesSpheres(1, 3);
	    Phase phase = new Phase();
	    
	    Potential2Group potential2 = new Potential2Group();
	    Potential2 p2 = new P2HardSphere(potential2);
	    potential2.setSpecies(speciesSpheres, speciesSpheres);
	    
	    Potential1 p1 = new P1TetheredHardSpheres();
	    p1.setSpecies(speciesSpheres);
	    
	    Controller controller = new Controller();
	    etomica.graphics.DisplayPhase displayPhase = new etomica.graphics.DisplayPhase();
        integratorHard.setTimeStep(0.01);
        
        //this method call invokes the mediator to tie together all the assembled components.
		Simulation.instance.elementCoordinator.go();
		                                    
        sim.makeAndDisplayFrame();
        
     //   controller.start();
    }//end of main
    
}//end of P1TetheredHardSpheres
   
