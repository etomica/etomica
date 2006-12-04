package etomica.modules.entropylottery;
import etomica.action.activity.ActivityIntegrate;
import etomica.integrator.IntegratorMC;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Space1D;
import etomica.species.SpeciesSpheresMono;

public class EntropyLottery extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono species;
    public Phase phase;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    
    public EntropyLottery(Space space) {
        super(space);
        defaults.makeLJDefaults();
        
        final int N = 30;  //number of atoms
        
        //controller and integrator
	    integrator = new IntegratorMC(this);
        MCMoveAtomAdjacent move = new MCMoveAtomAdjacent();
        integrator.getMoveManager().addMCMove(move);
        activityIntegrate = new ActivityIntegrate(this, integrator);
        activityIntegrate.setSleepPeriod(10);
        getController().addAction(activityIntegrate);

	    species = new SpeciesSpheresMono(this);
        
        //construct phase
	    phase = new Phase(this);
        phase.getAgent(species).setNMolecules(N);
        phase.setBoundary(new BoundaryRectangularNonperiodic(space));
        Vector dimensions = space.makeVector();
        dimensions.E(10);
        phase.getBoundary().setDimensions(dimensions);
//        new ConfigurationSequential(space).initializeCoordinates(phase);
        integrator.setPhase(phase);
		
//        PhaseImposePbc imposePBC = new PhaseImposePbc(phase);
//        integrator.addListener(new IntervalActionAdapter(imposePBC));
        
    }  
    
    public static void main(String[] args) {
        EntropyLottery sim = new EntropyLottery(Space1D.getInstance());
        sim.getController().actionPerformed();
    }
    
}


