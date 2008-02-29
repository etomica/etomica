package etomica.modules.entropylottery;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IVector;
import etomica.integrator.IntegratorMC;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space1d.Space1D;
import etomica.species.SpeciesSpheresMono;

public class EntropyLottery extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public SpeciesSpheresMono species;
    public Box box;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    
    public EntropyLottery(Space _space) {
        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        
        final int N = 30;  //number of atoms
        
        //controller and integrator
	    integrator = new IntegratorMC(this, potentialMaster);
        MCMoveAtomAdjacent move = new MCMoveAtomAdjacent(getRandom(), space);
        integrator.getMoveManager().addMCMove(move);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

	    species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        
        //construct box
	    box = new Box(new BoundaryRectangularNonperiodic(space, getRandom()), space);
        addBox(box);
        box.setNMolecules(species, N);
        IVector dimensions = space.makeVector();
        dimensions.E(10);
        box.getBoundary().setDimensions(dimensions);
        new ConfigurationZero(space).initializeCoordinates(box);
        integrator.setBox(box);
		
//        BoxImposePbc imposePBC = new BoxImposePbc(box);
//        integrator.addListener(new IntervalActionAdapter(imposePBC));
        
    }  
    
    public static void main(String[] args) {
        EntropyLottery sim = new EntropyLottery(Space1D.getInstance());
        sim.getController().actionPerformed();
    }
    
}


