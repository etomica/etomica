package etomica;

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
    
    public String getVersion() {return "P1IntraSimple:02.09.01/"+PotentialGroup.VERSION;}
    
    public Potential2 bonded;
    public Potential2 nonbonded;
    
    public P1IntraSimple() {
        this(Simulation.instance.hamiltonian.potential);
    }
    
    public P1IntraSimple(PotentialGroup parent) {
        super(parent);
    }
    
    /**
     * After constructing bonded potential with this group as its parent, this method
     * must be called with it as an argument to identify it as the bonded potential.
     */
    public void setBonded(Potential2 potential) {
        if(!this.contains(potential)) {
            throw new IllegalArgumentException("Error: Can identify only an existing child of P1IntraSimple as the bonded potential");
        }
        potential.setIterator(new ApiGeneral(parentSimulation().space,
	            new AtomIteratorList(),
	            new AtomIteratorBonds()));
    }
    
    /**
     * After constructing nonbonded potential with this group as its parent, this method
     * must be called with it as an argument to identify it as the nonbonded potential.
     */
    public void setNonbonded(Potential2 potential) {
        if(!this.contains(potential)) {
            throw new IllegalArgumentException("Error: Can identify only an existing child of P1IntraSimple as the nonbonded potential");
        }
        potential.setIterator(new ApiGeneral(parentSimulation().space,
	            new AtomIteratorList(),
	            new AtomIteratorNonbonded(parentSimulation())));
    }
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("General intramolecular potential with one bonded and one nonbonded potential");
        return info;
    }

    /**
     * Not implemented
     */
    public double energy(Atom a) {
        throw new RuntimeException("P1IntraSimple.energy method not implemented");
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
  //      Default.makeLJDefaults();//WHY DOESN'T THIS WORK?
        Default.TIME_STEP = etomica.units.LennardJones.Time.UNIT.fromSim(0.001);
        System.out.println("Time step: "+Default.TIME_STEP);
        Default.ATOM_SIZE = 1.0;
        Default.ATOM_MASS = 1.0;
        Default.POTENTIAL_WELL = 1.0;
        Default.TEMPERATURE = 1.0;//*/
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
	    
/*	    Potential2Group potential2 = new Potential2Group();
	    Potential2 p2 = new P2HardSphere(potential2);
	    potential2.setSpecies(speciesSpheres, speciesSpheres);
*/	    
	    P1IntraSimple p1 = new P1IntraSimple();
	    p1.setSpecies(speciesSpheres);
	    P2Fene p2Fene = new P2Fene(p1);
	    P2LennardJones p2LennardJones = new P2LennardJones(p1);
	    p1.setBonded(p2Fene);
	    p1.setNonbonded(p2LennardJones);
	    
	    Potential2Group p2 = new Potential2Group();
	    p2.setSpecies(speciesSpheres, speciesSpheres);
	    P2LennardJones p2LennardJones2 = new P2LennardJones(p2);
	    
	    Controller controller = new Controller();
	    etomica.graphics.DisplayPhase displayPhase = new etomica.graphics.DisplayPhase();
  //      integrator.setTimeStep(0.01);
  
        etomica.graphics.DeviceTrioControllerButton controlPanel = new etomica.graphics.DeviceTrioControllerButton();
        
        //this method call invokes the mediator to tie together all the assembled components.
		Simulation.instance.elementCoordinator.go();
		                                    
        sim.makeAndDisplayFrame();
        
     //   controller.start();
    }//end of main
    
}//end of P1TetheredHardSpheres
   
