package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeLeaf;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataAccumulator;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterCoupled;
import etomica.virial.ClusterWeight;
import etomica.virial.ConfigurationCluster;
import etomica.virial.MCMoveClusterAtomMulti;
import etomica.virial.MCMoveClusterMoleculeMulti;
import etomica.virial.MCMoveClusterReptateMulti;
import etomica.virial.MCMoveClusterRotateMoleculeMulti;
import etomica.virial.MCMoveClusterWiggleMulti;
import etomica.virial.MeterVirial;
import etomica.virial.P0Cluster;
import etomica.virial.BoxCluster;
import etomica.virial.SpeciesFactory;

/**
 * Generic simulation using Mayer sampling to evaluate cluster integrals
 */
public class SimulationVirial extends Simulation {


    /**
	 * Constructor for simulation to determine the ratio bewteen reference and target Clusters
	 */
	public SimulationVirial(Space space, SpeciesFactory speciesFactory, double temperature, ClusterWeight aSampleCluster, ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
		super(space,false);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        sampleCluster = aSampleCluster;
		int nMolecules = sampleCluster.pointCount();
		box = new BoxCluster(this,sampleCluster);
        box.getBoundary().setDimensions(Space.makeVector(new double[]{3.0,3.0,3.0}));
		species = speciesFactory.makeSpecies(this);
        getSpeciesManager().addSpecies(species);
        box.setNMolecules(species, nMolecules);
        
        if (refCluster instanceof ClusterCoupled) {
            ((ClusterCoupled)refCluster).setBox(box);
        }
        if (targetClusters[0] instanceof ClusterCoupled) {
            ((ClusterCoupled)targetClusters[0]).setBox(box);
        }

        integrator = new IntegratorMC(this, potentialMaster);
        // it's unclear what this accomplishes, but let's do it just for fun.
		integrator.setTemperature(temperature);
        integrator.setBox(box);
        integrator.getMoveManager().setEquilibrating(false);
        integrator.setEventInterval(1);
		ai = new ActivityIntegrate(integrator);
		getController().addAction(ai);
		
        if (species.getMoleculeType() instanceof AtomTypeLeaf) {
            mcMoveTranslate= new MCMoveClusterAtomMulti(this, potentialMaster, nMolecules-1);
        }
        else {
            mcMoveTranslate = new MCMoveClusterMoleculeMulti(potentialMaster,getRandom(),0.41,nMolecules-1);
            mcMoveRotate = new MCMoveClusterRotateMoleculeMulti(potentialMaster,getRandom(),nMolecules-1);
            mcMoveRotate.setStepSize(Math.PI);
            if (species instanceof SpeciesSpheres) {
                if (species.getMoleculeFactory().getNumChildAtoms() > 2) {
                    mcMoveWiggle = new MCMoveClusterWiggleMulti(this,potentialMaster, nMolecules);
                    integrator.getMoveManager().addMCMove(mcMoveWiggle);
                    mcMoveReptate = new MCMoveClusterReptateMulti(this,potentialMaster, nMolecules-1);
                    integrator.getMoveManager().addMCMove(mcMoveReptate);
                }
            }
            integrator.getMoveManager().addMCMove(mcMoveRotate);
        }
        integrator.getMoveManager().addMCMove(mcMoveTranslate);
		
		P0Cluster p0 = new P0Cluster(space);
		potentialMaster.addPotential(p0,new Species[]{});
		
        ConfigurationCluster configuration = new ConfigurationCluster(getRandom());
        configuration.initializeCoordinates(box);

        allValueClusters = new ClusterAbstract[targetClusters.length+1];
        allValueClusters[0] = refCluster;
        System.arraycopy(targetClusters,0,allValueClusters,1,targetClusters.length);
        setMeter(new MeterVirial(allValueClusters));
        ((MeterVirial)meter).setBox(box);
        setAccumulator(new AccumulatorRatioAverage());
	}
	
    private static final long serialVersionUID = 1L;
	public DataSource meter;
	public DataAccumulator accumulator;
	public DataPump accumulatorPump;
	public Species species;
	public ActivityIntegrate ai;
	public IntegratorMC integrator;
	public BoxCluster box;
    public ClusterAbstract[] allValueClusters;
    public ClusterWeight sampleCluster;
    public MCMoveBoxStep mcMoveTranslate;
    public MCMoveBoxStep mcMoveRotate;
    public MCMoveBoxStep mcMoveWiggle;
    public MCMoveBox mcMoveReptate;

	public void setMeter(DataSource newMeter) {
		meter = newMeter;
        if (accumulator != null) { 
            if (accumulatorPump != null) {
                integrator.removeIntervalAction(accumulatorPump);
                accumulatorPump = null;
            }
            setAccumulator(accumulator);
        }
	}

	public void setAccumulator(DataAccumulator newAccumulator) {
		accumulator = newAccumulator;
		if (accumulatorPump == null) {
			accumulatorPump = new DataPump(meter,accumulator);
            integrator.addIntervalAction(accumulatorPump);
		}
		else {
			accumulatorPump.setDataSink(accumulator);
		}
	}
}

