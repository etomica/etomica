package etomica.potential;

import etomica.EtomicaInfo;
import etomica.Simulation;
import etomica.Space;

/**
 * Intramolecular potential in which bonded and nonbonded atoms interact with a
 * hard potential (P2HardBond and P2HardSphere by default).
 * 
 * @author David Kofke
 */
 
public class P1BondedHardSpheres extends P1IntraSimple {

    public P1BondedHardSpheres() {
        this(Simulation.getDefault().space);
    }
    
    public P1BondedHardSpheres(Space space) {
        super(space, new P2HardBond(space), new P2HardSphere(space));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Bonded atoms tethered, nonbonded interact as hard spheres");
        return info;
    }

    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {
        etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic();
        Simulation.instance = sim;
	    IntegratorHard integratorHard = new IntegratorHard(sim.potentialMaster);
	    SpeciesSpheres speciesSpheres = new SpeciesSpheres(1, 3);
	    Phase phase = new Phase(sim.space);
	    
	    PotentialGroup potential2 = new PotentialGroup(2);
	    Potential2 p2 = new P2HardSphere();
	    sim.potentialMaster.setSpecies(p2, new Species[]{speciesSpheres});
	    //XXX need to add p2 to potential2 with an iterator
	    
	    PotentialGroup p1 = new P1BondedHardSpheres();
	    sim.potentialMaster.setSpecies(p1, new Species[]{speciesSpheres});
	    
	    Controller controller = new Controller();
	    etomica.graphics.DisplayPhase displayPhase = new etomica.graphics.DisplayPhase();
        integratorHard.setTimeStep(0.01);
        
        //this method call invokes the mediator to tie together all the assembled components.
		Simulation.instance.elementCoordinator.go();
		                                    
        sim.makeAndDisplayFrame();
        
     //   controller.start();
    }//end of main*/
    
}//end of P1TetheredHardSpheres
   
