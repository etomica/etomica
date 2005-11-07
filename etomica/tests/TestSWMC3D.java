package etomica.tests;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.config.ConfigurationFile;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.nbr.PotentialCalculationAgents;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.phase.Phase;
import etomica.potential.P2SquareWell;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
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
        super(space, false, new PotentialMasterCell(space, 1.5));
        defaults.makeLJDefaults();
        double sqwLambda = 1.5;
	    integrator = new IntegratorMC(this);
	    mcMoveAtom = new MCMoveAtom(this);
        mcMoveAtom.setStepSize(defaults.atomSize);
        integrator.addMCMove(mcMoveAtom);
        integrator.setEquilibrating(false);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setMaxSteps(500000);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(numAtoms);
	    phase = new Phase(this);
        phase.setDensity(0.7);
        integrator.addMCMoveListener(((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).makeMCMoveListener());
        potential = new P2SquareWell(space,defaults.atomSize,sqwLambda,defaults.potentialWell,false);
        ((PotentialMasterCell)potentialMaster).setNCells((int)(phase.getBoundary().dimensions().x(0)/potential.getRange()));
        ((PotentialMasterCell)potentialMaster).setRange(potential.getRange());
        potential.setCriterion(etomica.nbr.NeighborCriterion.ALL);
        potentialMaster.setSpecies(potential, new Species[] {species, species});
        
        new ConfigurationFile(space,"SWMC3D"+Integer.toString(numAtoms)).initializeCoordinates(phase);
        integrator.addPhase(phase);
        ((PotentialMasterCell)potentialMaster).calculate(phase, new PotentialCalculationAgents(potentialMaster));
        ((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).assignCellAll();
//        WriteConfiguration writeConfig = new WriteConfiguration("SWMC3D"+Integer.toString(numAtoms),phase,1);
//        integrator.addIntervalListener(writeConfig);
    }
 
    public static void main(String[] args) {
        int numAtoms = 500;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        TestSWMC3D sim = new TestSWMC3D(Space3D.getInstance(), numAtoms);
        sim.getDefaults().blockSize = 10;

        DataSource energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverage(sim);
        DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        new IntervalActionAdapter(energyManager, sim.integrator);
        
        sim.getController().actionPerformed();
        
        double avgPE = ((DataDoubleArray)((DataGroup)energyAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).getData()[0];
        avgPE /= numAtoms;
        System.out.println("PE/epsilon="+avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDoubleArray)((DataGroup)energyAccumulator.getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[0];
        Cv /= temp;
        Cv *= Cv/numAtoms;
        System.out.println("Cv/k="+Cv);
        
        if (Double.isNaN(avgPE) || Math.abs(avgPE+5.48) > 0.04) {
            System.exit(1);
        }
        // actual value ~0.56
        if (Double.isNaN(Cv) || Math.abs(Cv-0.7) > 0.6) {
            System.exit(1);
        }
    }
}