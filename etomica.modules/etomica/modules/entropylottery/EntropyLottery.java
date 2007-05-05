package etomica.modules.entropylottery;
import etomica.action.activity.ActivityIntegrate;
import etomica.integrator.IntegratorMC;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.IVector;
import etomica.space.Space;
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
        MCMoveAtomAdjacent move = new MCMoveAtomAdjacent(getRandom());
        integrator.getMoveManager().addMCMove(move);
        activityIntegrate = new ActivityIntegrate(this, integrator);
        getController().addAction(activityIntegrate);

	    species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        
        //construct phase
	    phase = new Phase(new BoundaryRectangularNonperiodic(space, getRandom()));
        addPhase(phase);
        phase.getAgent(species).setNMolecules(N);
        IVector dimensions = space.makeVector();
        dimensions.E(10);
        phase.getBoundary().setDimensions(dimensions);
        new ConfigurationZero().initializeCoordinates(phase);
        integrator.setPhase(phase);
		
//        PhaseImposePbc imposePBC = new PhaseImposePbc(phase);
//        integrator.addListener(new IntervalActionAdapter(imposePBC));
        
    }  
    
    public static void main(String[] args) {
        EntropyLottery sim = new EntropyLottery(Space1D.getInstance());
        sim.getController().actionPerformed();
    }
    
}


