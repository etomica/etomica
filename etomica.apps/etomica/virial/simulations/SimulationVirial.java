package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeLeaf;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataAccumulator;
import etomica.data.DataPump;
import etomica.data.meter.Meter;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.MCMove;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.Species;
import etomica.util.Default;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ConfigurationCluster;
import etomica.virial.IntegratorClusterMC;
import etomica.virial.MCMoveClusterAtom;
import etomica.virial.MCMoveClusterAtomMulti;
import etomica.virial.MCMoveClusterMolecule;
import etomica.virial.MCMoveClusterMoleculeMulti;
import etomica.virial.MCMoveClusterRotateMolecule3D;
import etomica.virial.MCMoveClusterRotateMoleculeMulti;
import etomica.virial.MeterVirial;
import etomica.virial.P0Cluster;
import etomica.virial.PhaseCluster;
import etomica.virial.SpeciesFactory;

/**
 * Generic simulation using Mayer sampling to evaluate cluster integrals
 */
public class SimulationVirial extends Simulation {

	/**
	 * Constructor for simulation to determine the ratio bewteen reference and target Clusters
	 */
	public SimulationVirial(Space space, Default defaults, SpeciesFactory speciesFactory, double temperature, ClusterWeight aSampleCluster, ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
		super(space,false,new PotentialMaster(space),Default.BIT_LENGTH,defaults);
        sampleCluster = aSampleCluster;
		int nMolecules = sampleCluster.pointCount();
		phase = new PhaseCluster(this,sampleCluster);
		species = speciesFactory.makeSpecies(this);//SpheresMono(this,AtomLinker.FACTORY);
        species.setNMolecules(nMolecules);
        phase.makeMolecules();
        
		integrator = new IntegratorClusterMC(this);
        // it's unclear what this accomplishes, but let's do it just for fun.
		integrator.setTemperature(temperature);
        integrator.setPhase(phase);
        integrator.setEquilibrating(false);
		ai = new ActivityIntegrate(this,integrator);
		ai.setInterval(1);
		getController().addAction(ai);
		
        if (species.getFactory().getType() instanceof AtomTypeLeaf) {
            if (false && (nMolecules == 2 || nMolecules > 5)) {
                System.out.println("using single-atom moves");
                mcMoveTranslate = new MCMoveClusterAtom(this);
                mcMoveTranslate.setStepSize(1.15);
            }
            else {
                System.out.println("using multi-atom moves");
                mcMoveTranslate= new MCMoveClusterAtomMulti(this, nMolecules-1);
                mcMoveTranslate.setStepSize(0.41);
            }
        }
        else {
            if (false && (nMolecules == 2 || nMolecules > 5)) {
                System.out.println("using single-atom moves");
                mcMoveTranslate = new MCMoveClusterMolecule(potentialMaster,3.0);
                mcMoveRotate = new MCMoveClusterRotateMolecule3D(potentialMaster,space);
                mcMoveRotate.setStepSize(Math.PI);
            }
            else {
                System.out.println("using multi-atom moves");
                mcMoveTranslate = new MCMoveClusterMoleculeMulti(potentialMaster,0.41,nMolecules-1);
                mcMoveRotate = new MCMoveClusterRotateMoleculeMulti(potentialMaster,space,nMolecules-1);
                mcMoveRotate.setStepSize(Math.PI);
            }
            integrator.getMoveManager().addMCMove(mcMoveRotate);
        }
        integrator.getMoveManager().addMCMove(mcMoveTranslate);
		
		P0Cluster p0 = new P0Cluster(space);
		potentialMaster.addPotential(p0,new Species[]{});
		
        ConfigurationCluster configuration = new ConfigurationCluster(space);
        configuration.setPhase(phase);
        configuration.initializeCoordinates(phase);

        allValueClusters = new ClusterAbstract[targetClusters.length+1];
        allValueClusters[0] = refCluster;
        System.arraycopy(targetClusters,0,allValueClusters,1,targetClusters.length);
        setMeter(new MeterVirial(allValueClusters,integrator));
        setAccumulator(new AccumulatorRatioAverage(this));
	}
	
    public IntervalActionAdapter accumulatorAA;
	public Meter meter;
	public DataAccumulator accumulator;
	public DataPump accumulatorPump;
	public Species species;
	public ActivityIntegrate ai;
	public IntegratorClusterMC integrator;
	public PhaseCluster phase;
    public ClusterAbstract[] allValueClusters;
    public ClusterWeight sampleCluster;
    public MCMove mcMoveTranslate;
    public MCMove mcMoveRotate;

	public void setMeter(Meter newMeter) {
		meter = newMeter;
		if (meter != null) {
			meter.setPhase(phase);
		}
        if (accumulator != null) { 
            if (accumulatorPump != null) {
                integrator.removeListener(accumulatorAA);
                accumulatorPump = null;
            }
            setAccumulator(accumulator);
        }
	}

	public void setAccumulator(DataAccumulator newAccumulator) {
		accumulator = newAccumulator;
		if (accumulatorPump == null) {
			accumulatorPump = new DataPump(meter,accumulator);
            accumulatorAA = new IntervalActionAdapter(accumulatorPump);
            integrator.addListener(accumulatorAA);
		}
		else {
			accumulatorPump.setDataSink(accumulator);
		}
	}
}

