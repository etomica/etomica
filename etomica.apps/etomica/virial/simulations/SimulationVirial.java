package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataAccumulator;
import etomica.data.DataPumpListener;
import etomica.data.IEtomicaDataSource;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;
import etomica.species.SpeciesSpheresRotating;
import etomica.virial.BoxCluster;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ConfigurationCluster;
import etomica.virial.MCMoveClusterAtomMulti;
import etomica.virial.MCMoveClusterMoleculeMulti;
import etomica.virial.MCMoveClusterReptateMulti;
import etomica.virial.MCMoveClusterRotateMoleculeMulti;
import etomica.virial.MCMoveClusterWiggleMulti;
import etomica.virial.MeterVirial;
import etomica.virial.P0Cluster;
import etomica.virial.SpeciesFactory;

/**
 * Generic simulation using Mayer sampling to evaluate cluster integrals
 */
public class SimulationVirial extends Simulation {


    /**
	 * Constructor for simulation to determine the ratio bewteen reference and target Clusters
	 */
	public SimulationVirial(ISpace space, SpeciesFactory speciesFactory, double temperature, ClusterWeight aSampleCluster, ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
	    this(space, speciesFactory, temperature, aSampleCluster, refCluster, targetClusters, false);
	}

	public SimulationVirial(ISpace space, SpeciesFactory speciesFactory, double temperature, ClusterWeight aSampleCluster, ClusterAbstract refCluster, ClusterAbstract[] targetClusters, boolean doWiggle) {
		super(space,false);
        PotentialMaster potentialMaster = new PotentialMaster();
        sampleCluster = aSampleCluster;
		int nMolecules = sampleCluster.pointCount();
		box = new BoxCluster(sampleCluster, space);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{3.0,3.0,3.0}));
		species = speciesFactory.makeSpecies(this, space);
        getSpeciesManager().addSpecies(species);
        box.setNMolecules(species, nMolecules);
        
        integrator = new IntegratorMC(this, potentialMaster);
        // it's unclear what this accomplishes, but let's do it just for fun.
		integrator.setTemperature(temperature);
        integrator.setBox(box);
        integrator.getMoveManager().setEquilibrating(false);
        integrator.setEventInterval(1);
		ai = new ActivityIntegrate(integrator);
		getController().addAction(ai);
		
		
        if (species instanceof SpeciesSpheresMono || species instanceof SpeciesSpheresRotating) {
            mcMoveTranslate = new MCMoveClusterAtomMulti(this, potentialMaster, space);
        }
        else {
            mcMoveTranslate = new MCMoveClusterMoleculeMulti(this, potentialMaster, space);
            mcMoveRotate = new MCMoveClusterRotateMoleculeMulti(potentialMaster,getRandom(), space);
            mcMoveRotate.setStepSize(Math.PI);
            if (species instanceof SpeciesSpheres) {
                if (doWiggle) {
                    mcMoveWiggle = new MCMoveClusterWiggleMulti(this,potentialMaster, nMolecules, space);
                    integrator.getMoveManager().addMCMove(mcMoveWiggle);
                    mcMoveReptate = new MCMoveClusterReptateMulti(this,potentialMaster, nMolecules-1);
                    integrator.getMoveManager().addMCMove(mcMoveReptate);
                }
            }
            integrator.getMoveManager().addMCMove(mcMoveRotate);
        }
        integrator.getMoveManager().addMCMove(mcMoveTranslate);
		
		P0Cluster p0 = new P0Cluster(space);
		potentialMaster.addPotential(p0,new ISpecies[]{});
		
        ConfigurationCluster configuration = new ConfigurationCluster(space);
        configuration.initializeCoordinates(box);

        allValueClusters = new ClusterAbstract[targetClusters.length+1];
        allValueClusters[0] = refCluster;
        System.arraycopy(targetClusters,0,allValueClusters,1,targetClusters.length);
        setMeter(new MeterVirial(allValueClusters));
        ((MeterVirial)meter).setBox(box);
        setAccumulator(new AccumulatorRatioAverage());
	}
	
    private static final long serialVersionUID = 1L;
	public IEtomicaDataSource meter;
	public DataAccumulator accumulator;
	public DataPumpListener accumulatorPump;
	public ISpecies species;
	public ActivityIntegrate ai;
	public IntegratorMC integrator;
	public BoxCluster box;
    public ClusterAbstract[] allValueClusters;
    public ClusterWeight sampleCluster;
    public MCMoveBoxStep mcMoveTranslate;
    public MCMoveBoxStep mcMoveRotate;
    public MCMoveBoxStep mcMoveWiggle;
    public MCMoveBox mcMoveReptate;

	public void setMeter(IEtomicaDataSource newMeter) {
		meter = newMeter;
        if (accumulator != null) { 
            if (accumulatorPump != null) {
                integrator.getEventManager().removeListener(accumulatorPump);
                accumulatorPump = null;
            }
            setAccumulator(accumulator);
        }
	}

	public void setAccumulator(DataAccumulator newAccumulator) {
		accumulator = newAccumulator;
		if (accumulatorPump == null) {
			accumulatorPump = new DataPumpListener(meter,accumulator);
            integrator.getEventManager().addListener(accumulatorPump);
		}
		else {
			accumulatorPump.setDataSink(accumulator);
		}
	}
}

