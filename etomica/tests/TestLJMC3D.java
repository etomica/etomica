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
import etomica.data.meter.MeterPressure;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.nbratom.cell.PotentialCalculationAgents;
import etomica.nbratom.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.space3d.Space3D;

/**
 * Simple Lennard-Jones Monte Carlo simulation in 3D.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 */
 
public class TestLJMC3D extends Simulation {
    
    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2LennardJones potential;
    public Controller controller;

    public TestLJMC3D(Space space, int numAtoms) {
        super(space,new PotentialMasterCell(space));
        space.setKinetic(false);
        Default.makeLJDefaults();
	    integrator = new IntegratorMC(potentialMaster);
	    mcMoveAtom = new MCMoveAtom(potentialMaster);
        mcMoveAtom.setStepSize(0.2*Default.ATOM_SIZE);
        integrator.addMCMove(mcMoveAtom);
        integrator.setEquilibrating(false);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(200000);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(numAtoms);
	    phase = new Phase(this);
        phase.setDensity(0.65);
        potential = new P2LennardJones(space);
        P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(potential);
        double truncationRadius = 3.0*potential.getSigma();
        if(truncationRadius > 0.5*phase.boundary().dimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is"+0.5*phase.boundary().dimensions().x(0));
        }
        potentialTruncated.setTruncationRadius(truncationRadius);
        ((PotentialMasterCell)potentialMaster).setNCells((int)(3.0*phase.boundary().dimensions().x(0)/potentialTruncated.getRange()));
        ((PotentialMasterCell)potentialMaster).setRange(potentialTruncated.getRange());
        potentialTruncated.setCriterion(etomica.nbr.NeighborCriterion.ALL);
        potentialMaster.setSpecies(potentialTruncated, new Species[] {species, species});
        integrator.addMCMoveListener(((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).makeMCMoveListener());
        
        ConfigurationFile config = new ConfigurationFile(space,"LJMC3D"+Integer.toString(numAtoms));
        phase.setConfiguration(config);
        integrator.addPhase(phase);
        ((PotentialMasterCell)potentialMaster).calculate(phase, new PotentialCalculationAgents());
        ((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).assignCellAll();
//        WriteConfiguration writeConfig = new WriteConfiguration("LJMC3D"+Integer.toString(numAtoms),phase,1);
//        integrator.addListener(writeConfig);
    }
 
    public static void main(String[] args) {
        int numAtoms = 500;
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        Default.BLOCK_SIZE = 10;
        TestLJMC3D sim = new TestLJMC3D(new Space3D(), numAtoms);

        MeterPressure pMeter = new MeterPressure(sim.potentialMaster,sim.space);
        pMeter.setTemperature(sim.integrator.getTemperature());
        pMeter.setPhase(sim.phase);
        AccumulatorAverage pAccumulator = new AccumulatorAverage();
        DataPump pPump = new DataPump(pMeter,pAccumulator);
        IntervalActionAdapter iaa = new IntervalActionAdapter(pPump,sim.integrator);
        iaa.setActionInterval(50);
        DataSource energyMeter = new IntegratorPotentialEnergy(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverage();
        DataPump energyManager = new DataPump(energyMeter,new DataSink[]{energyAccumulator});
        energyAccumulator.setBlockSize(50);
        new IntervalActionAdapter(energyManager, sim.integrator);
        
        sim.getController().actionPerformed();
        
        double[][] pData = (double[][])pAccumulator.getTranslator().fromArray(pAccumulator.getData()); 
        double Z = pData[AccumulatorAverage.AVERAGE.index][0]*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.getTemperature());
        double[] data = energyAccumulator.getData();
        double PE = data[AccumulatorAverage.AVERAGE.index]/numAtoms;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+PE);
        double temp = sim.integrator.getTemperature();
        double Cv = data[AccumulatorAverage.STANDARD_DEVIATION.index]/temp;
        Cv *= Cv/numAtoms;
        System.out.println("Cv/k="+Cv);
        
        if (Math.abs(Z-0.44) > 0.15) {
            System.exit(1);
        }
        if (Math.abs(PE+4.67) > 0.02) {
            System.exit(1);
        }
        // actual average closer ~ 0.61
        if (Math.abs(Cv-0.7) > 0.3) {
            System.exit(1);
        }
    }

}
