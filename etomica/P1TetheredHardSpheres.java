package etomica;

/**
 * Intramolecular potential in which bonded atoms interact with a hard tether
 * potential, and nonbonded atoms interact as hard spheres.
 *
 * @author David Kofke
 */
 
public class P1TetheredHardSpheres extends PotentialGroup implements Potential1.Intramolecular {
    
    public final P2HardSphere p2HardSphere;
    public final P2Tether p2Tether;
    
    public P1TetheredHardSpheres() {
        this(Simulation.getDefault().space);
    }
    
    public P1TetheredHardSpheres(Space space) {
        super(1, space);
        p2HardSphere = new P2HardSphere();
        p2Tether = new P2Tether();
        addPotential(p2Tether, new ApiInnerVariable(new AtomIteratorList(),
	            new AtomIteratorBonds()));
	    addPotential(p2HardSphere, new ApiInnerVariable(new AtomIteratorList(),
	            new AtomIteratorNonbonded(parent.simulation())));
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
	    
	    PotentialGroup potential2 = new PotentialGroup(2);
	    Potential2 p2 = new P2HardSphere();
	    sim.potentialMaster.setSpecies(p2, new Species[]{speciesSpheres});
	    //XXX need to add p2 to potential2 with an iterator
	    
	    PotentialGroup p1 = new P1TetheredHardSpheres();
	    sim.potentialMaster.setSpecies(p1, new Species[]{speciesSpheres});
	    
	    Controller controller = new Controller();
	    etomica.graphics.DisplayPhase displayPhase = new etomica.graphics.DisplayPhase();
        integratorHard.setTimeStep(0.01);
        
        //this method call invokes the mediator to tie together all the assembled components.
		Simulation.instance.elementCoordinator.go();
		                                    
        sim.makeAndDisplayFrame();
        
     //   controller.start();
    }//end of main
    
}//end of P1TetheredHardSpheres
   
