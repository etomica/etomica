package etomica.potential;

import etomica.EtomicaInfo;
import etomica.Simulation;
import etomica.Space;
import etomica.atom.iterator.ApiBuilder;


/**
 * Generic intramolecular potential group, having one potential for bonded
 * atoms, and a different potential for unbonded ones.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 09/01/02 (DAK) new
  */
 
public class P1IntraSimple extends PotentialGroup implements Potential1.Intramolecular {
    
    public Potential2 bonded;
    public Potential2 nonBonded;
    
    public P1IntraSimple() {
        this(Simulation.getDefault().space);
    }
    
    public P1IntraSimple(Space space) {
        super(1, space);
    }
    
    public P1IntraSimple(Space space, Potential2 bonded, Potential2 nonbonded) {
        this(space);
        setBonded(bonded);
        setNonbonded(nonbonded);
    }
    
    /**
     * After constructing bonded potential with this group as its parent, this method
     * must be called with it as an argument to identify it as the bonded potential.
     */
    public void setBonded(Potential2 potential) {
        if(potential == null) return;
        if(bonded != null) removePotential(bonded);
        bonded = potential;
        addPotential(potential, ApiBuilder.makeAdjacentPairIterator());
    }
    
    /**
     * After constructing nonbonded potential with this group as its parent, this method
     * must be called with it as an argument to identify it as the nonbonded potential.
     */
    public void setNonbonded(Potential2 potential) {
        if(potential == null) return;
        if(nonBonded != null) removePotential(nonBonded);
        nonBonded = potential;
        addPotential(potential, ApiBuilder.makeNonAdjacentPairIterator());
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("General intramolecular potential with one bonded and one nonbonded potential");
        return info;
    }

    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {
  //      Default.makeLJDefaults();//WHY DOESN'T THIS WORK?
        Default.TIME_STEP = etomica.units.systems.LJ.SYSTEM.time().toSim(0.001);
        System.out.println("Time step: "+Default.TIME_STEP);
        Default.ATOM_SIZE = 1.0;
        Default.ATOM_MASS = 1.0;
        Default.POTENTIAL_WELL = 1.0;
        Default.TEMPERATURE = 1.0;
        System.out.println(Default.ATOM_SIZE);
        System.out.println(Default.ATOM_MASS);
        System.out.println(Default.POTENTIAL_WELL);
        System.out.println(Default.TEMPERATURE);
//        etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new Space3D());
        etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic();
        etomica.graphics.DefaultGraphic.ATOM_COLOR = java.awt.Color.red;
        Default.TRUNCATE_POTENTIALS = false;
        Simulation.instance = sim;
	    IntegratorVelocityVerlet integrator = new IntegratorVelocityVerlet();
//	    integrator.setIsothermal(true);
	    integrator.setSleepPeriod(1);
	    SpeciesSpheres speciesSpheres = new SpeciesSpheres(15, 5);
	    Phase phase = new Phase();
	    phase.setLrcEnabled(false);
	    
	    Potential2Group potential2 = new Potential2Group();
	    Potential2 p2 = new P2HardSphere(potential2);
	    potential2.setSpecies(speciesSpheres, speciesSpheres);

//	    P1IntraSimple p1 = new P1IntraSimple();
//	    p1.setSpecies(new Species[] {speciesSpheres});
//	    P2Fene p2Fene = new P2Fene(p1);
//	    P2LennardJones p2LennardJones = new P2LennardJones(p1);
//	    p1.setBonded(p2Fene);
//	    p1.setNonbonded(p2LennardJones);
	    
	    PotentialGroup p2 = new PotentialGroup(2);
	    
	    sim.potentialMaster.setSpecies(p2, new Species[] {speciesSpheres});
	    P2LennardJones p2LennardJones2 = new P2LennardJones();
	    //XXX p2LennardJones2 needs to be added to p2 with an iterator!
	    
	    Controller controller = new Controller();
	    etomica.graphics.DisplayPhase displayPhase = new etomica.graphics.DisplayPhase();
  //      integrator.setTimeStep(0.01);
  
        etomica.graphics.DeviceTrioControllerButton controlPanel = new etomica.graphics.DeviceTrioControllerButton();
        
        //this method call invokes the mediator to tie together all the assembled components.
		Simulation.instance.elementCoordinator.go();
		                                    
        sim.makeAndDisplayFrame();
        
     //   controller.start();
    }//end of main*/
    
}//end of P1TetheredHardSpheres
   
