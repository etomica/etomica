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
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.nbratom.cell.PotentialCalculationAgents;
import etomica.nbratom.cell.PotentialMasterCell;
import etomica.potential.P2SquareWell;
import etomica.space3d.Space3D;

/**
 * Simple hard-sphere Monte Carlo simulation in 2D.
 *
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
	    integrator = new IntegratorMC(potentialMaster);
	    mcMoveAtom = new MCMoveAtom(potentialMaster);
        mcMoveAtom.setStepSize(2*Default.ATOM_SIZE);
        integrator.addMCMove(mcMoveAtom);
        integrator.setEquilibrating(false);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(100000000/numAtoms);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(numAtoms);
	    phase = new Phase(this);
        phase.setDensity(0.1);
        integrator.addMCMoveListener(((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).makeMCMoveListener());
        potential = new P2SquareWell(space);
        ((PotentialMasterCell)potentialMaster).setNCells((int)(phase.boundary().dimensions().x(0)/potential.getRange()));
        potential.setCriterion(etomica.nbr.NeighborCriterion.ALL);
        potentialMaster.setSpecies(potential, new Species[] {species, species});
        
        ConfigurationFile config = new ConfigurationFile(space,"LJMC3D"+Integer.toString(numAtoms));
//        phase.setConfiguration(config);
        integrator.addPhase(phase);
        ((PotentialMasterCell)potentialMaster).calculate(phase, new PotentialCalculationAgents());
        ((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).assignCellAll();
//        integrator.addIntervalListener(new PhaseImposePbc(phase));
//        WriteConfiguration writeConfig = new WriteConfiguration("LJMC3D"+Integer.toString(numAtoms),phase,1);
//        integrator.addIntervalListener(writeConfig);
    }
 
    public static void main(String[] args) {
        int numAtoms = 500;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        Default.BLOCK_SIZE = 10;
        TestSWMC3D sim = new TestSWMC3D(new Space3D(), numAtoms);

        SimulationGraphic graphic = new SimulationGraphic(sim);
        graphic.makeAndDisplayFrame();
        ColorSchemeByType.setColor(sim.species.getFactory().getType(), java.awt.Color.red);
        
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
        
        if (Math.abs(PE+4.52) > 0.03) {
            System.exit(1);
        }
        if (Math.abs(Cv-0.034) > 0.01) {
            System.exit(1);
        }
    }

}
