package etomica.modules.ljmd;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationSequential;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.phase.Phase;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class Ljmd extends Simulation {
    
    public SpeciesSpheresMono species;
    public Phase phase;
    public IntegratorVelocityVerlet integrator;
    public ActivityIntegrate activityIntegrate;
    
    public Ljmd(Space space) {
        super(space);
        defaults.makeLJDefaults();
        defaults.boxSize = 14.0;
        
        final int N = 182;  //number of atoms
        
        
        //controller and integrator
	    integrator = new IntegratorVelocityVerlet(this);
	    integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        activityIntegrate = new ActivityIntegrate(this, integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        integrator.setTimeStep(0.01);
     //   integrator.setDoSleep(false);

	    //species and potentials
	//    SpeciesSphereWells disks = new SpeciesSphereWells(this);//index 1
	    species = new SpeciesSpheresMono(this);//index 1
	    species.setName("");
	    species.setNMolecules(N);
        
        //instantiate several potentials for selection in combo-box
	    P2LennardJones potential = new P2LennardJones(this);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(potential,2.5);
	    potentialMaster.addPotential(p2Truncated, new Species[]{species, species});
	    
        //construct phase
	    phase = new Phase(this);
        new ConfigurationSequential(space).initializeCoordinates(phase);
        integrator.setPhase(phase);
		
        PhaseImposePbc imposePBC = new PhaseImposePbc(phase);
        integrator.addListener(new IntervalActionAdapter(imposePBC));
        
    }//end of constructor    
    
    public static void main(String[] args) {
        Space space = Space2D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
            
        Ljmd sim = new Ljmd(space);
        sim.getController().actionPerformed();
    }//end of main
    
}


