package etomica.tests;
import etomica.ConfigurationFile;
import etomica.Controller;
import etomica.DataSink;
import etomica.DataSource;
import etomica.Default;
import etomica.IntegratorPotentialEnergy;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.nbratom.cell.PotentialCalculationAgents;
import etomica.nbratom.cell.PotentialMasterCell;
import etomica.potential.P2SquareWell;
import etomica.space3d.Space3D;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D.
 * Initial configurations at http://gordon.eng.buffalo.edu/etomica/tests/
 * @author David Kofke
 */
 
public class TestSWMC3D extends Simulation {
    
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2SquareWell potential;
    public Controller controller;

    public TestSWMC3D(Space space, int numAtoms) {
        super(space, new PotentialMasterCell(space));
        space.setKinetic(false);
        Default.makeLJDefaults();
        double sqwLambda = 1.5;
	    integrator = new IntegratorMC(potentialMaster);
	    mcMoveAtom = new MCMoveAtom(potentialMaster);
        mcMoveAtom.setStepSize(Default.ATOM_SIZE);
        integrator.addMCMove(mcMoveAtom);
        integrator.setEquilibrating(false);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(500000);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(numAtoms);
	    phase = new Phase(this);
        phase.setDensity(0.7);
        integrator.addMCMoveListener(((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).makeMCMoveListener());
        potential = new P2SquareWell(space,Default.ATOM_SIZE,sqwLambda,Default.POTENTIAL_WELL);
        ((PotentialMasterCell)potentialMaster).setNCells((int)(phase.boundary().dimensions().x(0)/potential.getRange()));
        ((PotentialMasterCell)potentialMaster).setRange(potential.getRange());
        potential.setCriterion(etomica.nbr.NeighborCriterion.ALL);
        potentialMaster.setSpecies(potential, new Species[] {species, species});
        
        ConfigurationFile config = new ConfigurationFile(space,"SWMC3D"+Integer.toString(numAtoms));
        phase.setConfiguration(config);
        integrator.addPhase(phase);
        ((PotentialMasterCell)potentialMaster).calculate(phase, new PotentialCalculationAgents());
        ((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).assignCellAll();
//        WriteConfiguration writeConfig = new WriteConfiguration("SWMC3D"+Integer.toString(numAtoms),phase,1);
//        integrator.addIntervalListener(writeConfig);
    }
 
    public static void main(String[] args) {
        int numAtoms = 500;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        Default.BLOCK_SIZE = 10;
        TestSWMC3D sim = new TestSWMC3D(new Space3D(), numAtoms);

        DataSource energyMeter = new IntegratorPotentialEnergy(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverage();
        DataPump energyManager = new DataPump(energyMeter,new DataSink[]{energyAccumulator});
        energyAccumulator.setBlockSize(50);
        new IntervalActionAdapter(energyManager, sim.integrator);
        
        sim.getController().actionPerformed();
        
        double[] data = energyAccumulator.getData();
        double PE = data[AccumulatorAverage.AVERAGE.index]/numAtoms;
        System.out.println("PE/epsilon="+PE);
        double temp = sim.integrator.temperature();
        double Cv = data[AccumulatorAverage.STANDARD_DEVIATION.index]/(temp*temp*numAtoms);
        System.out.println("Cv/k="+Cv);
        
        if (Math.abs(PE+5.48) > 0.04) {
            System.exit(1);
        }
        if (Math.abs(Cv-0.031) > 0.012) {
            System.exit(1);
        }
    }
}