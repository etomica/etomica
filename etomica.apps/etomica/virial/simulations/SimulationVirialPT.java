package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataAccumulator;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.DataSourceAcceptanceProbability;
import etomica.data.DataSourceAcceptanceRatio;
import etomica.data.meter.Meter;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorPT;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMovePhaseStep;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.Species;
import etomica.util.Default;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ConfigurationCluster;
import etomica.virial.MCMoveClusterAtom;
import etomica.virial.MCMoveClusterAtomMulti;
import etomica.virial.MCMoveClusterMolecule;
import etomica.virial.MCMoveClusterMoleculeMulti;
import etomica.virial.MCMoveClusterRotateMolecule3D;
import etomica.virial.MeterVirial;
import etomica.virial.P0Cluster;
import etomica.virial.PhaseCluster;
import etomica.virial.SpeciesFactory;
import etomica.virial.paralleltempering.MCMoveSwapCluster;

/**
 * Generic simulation using Mayer sampling to evaluate cluster integrals
 */
public class SimulationVirialPT extends Simulation {

	/**
	 * Constructor for simulation to determine the ratio bewteen reference and target Clusters
	 */
	public SimulationVirialPT(Space space, Default defaults, SpeciesFactory speciesFactory, 
			double[] temperature, ClusterWeight.Factory sampleClusterFactory, 
			ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
		super(space,false,new PotentialMaster(space),Default.BIT_LENGTH,defaults);
        
		int nMolecules = refCluster.pointCount();
		species = speciesFactory.makeSpecies(this);//SpheresMono(this,AtomLinker.FACTORY);
        getSpeciesRoot().addSpecies(species);

        phase = new PhaseCluster[temperature.length];
        integrator = new IntegratorMC[temperature.length];
        meter = new Meter[temperature.length];
        accumulator = new DataAccumulator[temperature.length];
        accumulatorPump = new DataPump[temperature.length];
        dumb = new IntervalActionAdapter[temperature.length];
        mcMoveAtom1 = new MCMoveAtom[temperature.length];
        mcMoveMulti = new MCMovePhaseStep[temperature.length];
        mcMoveRotate = new MCMovePhaseStep[temperature.length];
        meterAccept = new DataSource[temperature.length-1];
        meterAcceptP = new DataSource[temperature.length-1];
        
        // Parallel tempering would sorta work without separate instances of the clusters
        // but the value caching based on coordinate pair set ID would be confused because
        // the coordinate pair sets from different phases could have the same ID and different
        // values
        sampleCluster = new ClusterWeight[temperature.length];
        allValueClusters = new ClusterAbstract[temperature.length][targetClusters.length+1];
        integratorPT = new IntegratorPT(potentialMaster,getRandom(),MCMoveSwapCluster.FACTORY);
//        integratorPT.setSwapInterval(2);
        ai = new ActivityIntegrate(this,integratorPT);
        ai.setInterval(1);
        getController().addAction(ai);
        
        for (int iTemp=0; iTemp<temperature.length; iTemp++) {
            allValueClusters[iTemp][0] = refCluster.makeCopy();
            allValueClusters[iTemp][0].setTemperature(temperature[iTemp]);
            for (int i=0; i<targetClusters.length; i++) {
                allValueClusters[iTemp][i+1] = targetClusters[i].makeCopy();
                allValueClusters[iTemp][i+1].setTemperature(temperature[iTemp]);
            }

            sampleCluster[iTemp] = sampleClusterFactory.makeWeightCluster(allValueClusters[iTemp]);
            sampleCluster[iTemp].setTemperature(temperature[iTemp]);
            phase[iTemp] = new PhaseCluster(this,sampleCluster[iTemp]);
            phase[iTemp].getAgent(species).setNMolecules(nMolecules);
            
            integrator[iTemp] = new IntegratorMC(this);
            integrator[iTemp].setTemperature(temperature[iTemp]);
            integrator[iTemp].setPhase(phase[iTemp]);
            integrator[iTemp].setEquilibrating(false);
            integratorPT.addIntegrator(integrator[iTemp]);
            
            MCMoveManager moveManager = integrator[iTemp].getMoveManager();
            
            if (phase[iTemp].molecule(0).isLeaf()) {
                mcMoveAtom1[iTemp] = new MCMoveClusterAtom(this);
                mcMoveAtom1[iTemp].setStepSize(1.15);
                moveManager.addMCMove(mcMoveAtom1[iTemp]);
                if (nMolecules>2) {
                    mcMoveMulti[iTemp] = new MCMoveClusterAtomMulti(this, nMolecules-1);
                    mcMoveMulti[iTemp].setStepSize(0.41);
                    moveManager.addMCMove(mcMoveMulti[iTemp]);
                }
            }
            else {
                mcMoveAtom1[iTemp] = new MCMoveClusterMolecule(potentialMaster, getRandom(), 3.0);
                moveManager.addMCMove(mcMoveAtom1[iTemp]);
                mcMoveRotate[iTemp] = new MCMoveClusterRotateMolecule3D(potentialMaster,getRandom());
                mcMoveRotate[iTemp].setStepSize(Math.PI);
                moveManager.addMCMove(mcMoveRotate[iTemp]);
                if (nMolecules>2) {
                    mcMoveMulti[iTemp] = new MCMoveClusterMoleculeMulti(potentialMaster, getRandom(), 0.41, nMolecules-1);
                    moveManager.addMCMove(mcMoveMulti[iTemp]);
                }
            }
            
            ConfigurationCluster configuration = new ConfigurationCluster(getRandom());
            configuration.initializeCoordinates(phase[iTemp]);
            
            setMeter(iTemp,new MeterVirial(allValueClusters[iTemp]));
//            meter[iTemp].getDataInfo().setLabel("Target/Refernce Ratio "+iTemp);
            setAccumulator(iTemp,new AccumulatorRatioAverage(getDefaults().blockSize));

            if(iTemp>0) {
                MCMove swapMove = integratorPT.getMoveManager().getMCMoves()[iTemp-1];
                meterAccept[iTemp-1] = new DataSourceAcceptanceRatio(swapMove);
                meterAcceptP[iTemp-1] = new DataSourceAcceptanceProbability(swapMove);
            }
        }
        P0Cluster p0 = new P0Cluster(space);
        potentialMaster.addPotential(p0,new Species[]{});
	}
	
    private static final long serialVersionUID = 1L;
	public DataSource[] meter;
    public DataSource[] meterAccept;
    public DataSource[] meterAcceptP;
	public DataAccumulator[] accumulator;
	public DataPump[] accumulatorPump;
    public IntervalActionAdapter[] dumb;
	public Species species;
	public ActivityIntegrate ai;
	public IntegratorMC[] integrator;
	public PhaseCluster[] phase;
    public ClusterAbstract[][] allValueClusters;
    public ClusterWeight[] sampleCluster;
    public MCMoveAtom[] mcMoveAtom1;
    public MCMovePhaseStep[] mcMoveRotate;
    public MCMovePhaseStep[] mcMoveMulti;
    public IntegratorPT integratorPT;

	public void setMeter(int i, DataSource newMeter) {
		if (accumulator[i] != null) { 
			if (accumulatorPump[i] != null) {
                integrator[i].removeListener(dumb[i]);
				accumulatorPump[i] = null;
                dumb[i] = null;
			}
			accumulator[i] = null;
		}
		meter[i] = newMeter;
		if (meter[i] != null && meter[i] instanceof Meter) {
			((Meter)meter[i]).setPhase(phase[i]);
		}
	}

	public void setAccumulator(int i, DataAccumulator newAccumulator) {
		accumulator[i] = newAccumulator;
		if (accumulatorPump[i] == null) {
			accumulatorPump[i] = new DataPump(meter[i],accumulator[i]);
		}
		else {
			accumulatorPump[i].setDataSink(accumulator[i]);
		}
        dumb[i] = new IntervalActionAdapter(accumulatorPump[i]);
		integrator[i].addListener(dumb[i]);
	}
}

