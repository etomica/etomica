package etomica.virial.simulations;

import etomica.Simulation;
import etomica.Species;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataAccumulator;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.DataSourceAcceptanceProbability;
import etomica.data.DataSourceAcceptanceRatio;
import etomica.data.meter.Meter;
import etomica.integrator.IntegratorPT;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.MCMove;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ConfigurationCluster;
import etomica.virial.IntegratorClusterMC;
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
	public SimulationVirialPT(Space space, SpeciesFactory speciesFactory, double[] temperature, ClusterWeight.Factory sampleClusterFactory, ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
		super(space);
        
		int nMolecules = refCluster.pointCount();
		species = speciesFactory.makeSpecies(this);//SpheresMono(this,AtomLinker.FACTORY);
        species.setNMolecules(nMolecules);

        phase = new PhaseCluster[temperature.length];
        integrator = new IntegratorClusterMC[temperature.length];
        meter = new Meter[temperature.length];
        accumulator = new DataAccumulator[temperature.length];
        accumulatorPump = new DataPump[temperature.length];
        dumb = new IntervalActionAdapter[temperature.length];
        mcMoveAtom1 = new MCMoveAtom[temperature.length];
        mcMoveMulti = new MCMove[temperature.length];
        mcMoveRotate = new MCMove[temperature.length];
        meterAccept = new DataSource[temperature.length-1];
        meterAcceptP = new DataSource[temperature.length-1];
        
        // Parallel tempering would sorta work without separate instances of the clusters
        // but the value caching based on coordinate pair set ID would be confused because
        // the coordinate pair sets from different phases could have the same ID and different
        // values
        sampleCluster = new ClusterWeight[temperature.length];
        allValueClusters = new ClusterAbstract[temperature.length][targetClusters.length+1];
        integratorPT = new IntegratorPT(potentialMaster,MCMoveSwapCluster.FACTORY);
//        integratorPT.setSwapInterval(2);
        ai = new ActivityIntegrate(integratorPT);
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
            phase[iTemp].makeMolecules();
            
            integrator[iTemp] = new IntegratorClusterMC(potentialMaster);
            integrator[iTemp].setTemperature(temperature[iTemp]);
            integrator[iTemp].addPhase(phase[iTemp]);
            integrator[iTemp].setEquilibrating(false);
            integratorPT.addIntegrator(integrator[iTemp]);
            
            if (phase[iTemp].randomMolecule().node.isLeaf()) {
                mcMoveAtom1[iTemp] = new MCMoveClusterAtom(potentialMaster);
                mcMoveAtom1[iTemp].setStepSize(1.15);
                integrator[iTemp].addMCMove(mcMoveAtom1[iTemp]);
                if (nMolecules>2) {
                    mcMoveMulti[iTemp] = new MCMoveClusterAtomMulti(potentialMaster, nMolecules-1);
                    mcMoveMulti[iTemp].setStepSize(0.41);
                    integrator[iTemp].addMCMove(mcMoveMulti[iTemp]);
                }
            }
            else {
                mcMoveAtom1[iTemp] = new MCMoveClusterMolecule(potentialMaster);
                mcMoveAtom1[iTemp].setStepSize(3.0);
                integrator[iTemp].addMCMove(mcMoveAtom1[iTemp]);
                mcMoveRotate[iTemp] = new MCMoveClusterRotateMolecule3D(potentialMaster,space);
                mcMoveRotate[iTemp].setStepSize(Math.PI);
                integrator[iTemp].addMCMove(mcMoveRotate[iTemp]);
                if (nMolecules>2) {
                    mcMoveMulti[iTemp] = new MCMoveClusterMoleculeMulti(potentialMaster, nMolecules-1);
                    mcMoveMulti[iTemp].setStepSize(0.41);
                    integrator[iTemp].addMCMove(mcMoveMulti[iTemp]);
                }
            }
            
            ConfigurationCluster configuration = new ConfigurationCluster(space);
            configuration.setPhase(phase[iTemp]);
            configuration.initializeCoordinates(phase[iTemp]);
            
            setMeter(iTemp,new MeterVirial(allValueClusters[iTemp],integrator[iTemp]));
//            meter[iTemp].getDataInfo().setLabel("Target/Refernce Ratio "+iTemp);
            setAccumulator(iTemp,new AccumulatorRatioAverage());

            if(iTemp>0) {
                meterAccept[iTemp-1] = new DataSourceAcceptanceRatio(integratorPT.swapMoves()[iTemp-1]);
                meterAcceptP[iTemp-1] = new DataSourceAcceptanceProbability(integratorPT.swapMoves()[iTemp-1]);
            }
        }
        P0Cluster p0 = new P0Cluster(space);
        potentialMaster.setSpecies(p0,new Species[]{});
	}
	
	public Meter[] meter;
    public DataSource[] meterAccept;
    public DataSource[] meterAcceptP;
	public DataAccumulator[] accumulator;
	public DataPump[] accumulatorPump;
    public IntervalActionAdapter[] dumb;
	public Species species;
	public ActivityIntegrate ai;
	public IntegratorClusterMC[] integrator;
	public PhaseCluster[] phase;
    public ClusterAbstract[][] allValueClusters;
    public ClusterWeight[] sampleCluster;
    public MCMoveAtom[] mcMoveAtom1;
    public MCMove[] mcMoveRotate;
    public MCMove[] mcMoveMulti;
    public IntegratorPT integratorPT;

	public void setMeter(int i, Meter newMeter) {
		if (accumulator[i] != null) { 
			if (accumulatorPump[i] != null) {
                integrator[i].removeListener(dumb[i]);
				accumulatorPump[i] = null;
                dumb[i] = null;
			}
			accumulator[i] = null;
		}
		meter[i] = newMeter;
		if (meter[i] != null) {
			meter[i].setPhase(phase[i]);
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

