package etomica;

/**
 * Intramolecular potential in which bonded and nonbonded atoms interact with a
 * hard potential (P2HardBond and P2HardSphere by default).
 * 
 * @author David Kofke
 */
 
public class P1BondedHardSpheres extends PotentialGroup implements Potential1.Intramolecular {
    
    private Potential2 p2HardNonbonded;
    private Potential2 p2HardBond;
    
    public P1BondedHardSpheres() {
        this(Simulation.getDefault().space);
    }
    
    public P1BondedHardSpheres(Space space) {
        super(1, space);
        p2HardBond= new P2HardBond(space);
        addPotential(p2HardBond, new ApiIntragroup(new ApiInnerVariable(new AtomIteratorBasis(),
                new AtomIteratorSequencerBonded())));
        p2HardNonbonded = new P2HardSphere(space);
        ApiIntragroup nonBonded = new ApiIntragroup();
        nonBonded.getInnerIterator().setNumToSkip(2);
	    addPotential(p2HardNonbonded, nonBonded);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Bonded atoms tethered, nonbonded interact as hard spheres");
        return info;
    }
    
    public void setBonded(Potential2 potential) {
        removePotential(p2HardBond);
        p2HardBond = potential;
        addPotential(potential, new ApiIntragroup(new ApiInnerVariable(new AtomIteratorBasis(),
                new AtomIteratorSequencerBonded())));
    }
    
    public void setNonBondedPotential(Potential2 potential) {
        removePotential(p2HardNonbonded);
        p2HardNonbonded = potential;
        ApiIntragroup nonBonded = new ApiIntragroup();
        nonBonded.getInnerIterator().setNumToSkip(2);
        addPotential(potential, nonBonded);
    }
    
    public double energy(Atom a) {
        return 0.0;
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
   
