package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtomGroup;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataAccumulator;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.data.DataSourceAcceptanceProbability;
import etomica.data.DataSourceAcceptanceRatio;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorPT;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.Species;
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
import etomica.virial.BoxCluster;
import etomica.virial.SpeciesFactory;
import etomica.virial.paralleltempering.MCMoveSwapCluster;

/**
 * Generic simulation using Mayer sampling to evaluate cluster integrals
 */
public class SimulationVirialPT extends Simulation {

	/**
	 * Constructor for simulation to determine the ratio bewteen reference and target Clusters
	 */
	public SimulationVirialPT(Space space, SpeciesFactory speciesFactory, 
			double[] temperature, ClusterWeight.Factory sampleClusterFactory, 
			ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
		super(space,false);
        PotentialMaster potentialMaster = new PotentialMaster(space);
		int nMolecules = refCluster.pointCount();
		species = speciesFactory.makeSpecies(this);//SpheresMono(this,AtomLinker.FACTORY);
        getSpeciesManager().addSpecies(species);

        box = new BoxCluster[temperature.length];
        integrator = new IntegratorMC[temperature.length];
        meter = new MeterVirial[temperature.length];
        accumulator = new DataAccumulator[temperature.length];
        accumulatorPump = new DataPump[temperature.length];
        mcMoveAtom1 = new MCMoveAtom[temperature.length];
        mcMoveMulti = new MCMoveBoxStep[temperature.length];
        mcMoveRotate = new MCMoveBoxStep[temperature.length];
        meterAccept = new DataSource[temperature.length-1];
        meterAcceptP = new DataSource[temperature.length-1];
        
        // Parallel tempering would sorta work without separate instances of the clusters
        // but the value caching based on coordinate pair set ID would be confused because
        // the coordinate pair sets from different boxs could have the same ID and different
        // values
        sampleCluster = new ClusterWeight[temperature.length];
        allValueClusters = new ClusterAbstract[temperature.length][targetClusters.length+1];
        integratorPT = new IntegratorPT(getRandom(),MCMoveSwapCluster.FACTORY);
//        integratorPT.setSwapInterval(2);
        integratorPT.setEventInterval(1);
        ai = new ActivityIntegrate(integratorPT);
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
            box[iTemp] = new BoxCluster(this,sampleCluster[iTemp]);
            box[iTemp].setNMolecules(species, nMolecules);
            
            integrator[iTemp] = new IntegratorMC(this, potentialMaster);
            integrator[iTemp].setTemperature(temperature[iTemp]);
            integrator[iTemp].setBox(box[iTemp]);
            integrator[iTemp].getMoveManager().setEquilibrating(false);
            integratorPT.addIntegrator(integrator[iTemp]);
            
            MCMoveManager moveManager = integrator[iTemp].getMoveManager();
            
            if (!(box[iTemp].molecule(0) instanceof IAtomGroup)) {
                mcMoveAtom1[iTemp] = new MCMoveClusterAtom(this, potentialMaster);
                mcMoveAtom1[iTemp].setStepSize(1.15);
                moveManager.addMCMove(mcMoveAtom1[iTemp]);
                if (nMolecules>2) {
                    mcMoveMulti[iTemp] = new MCMoveClusterAtomMulti(this, potentialMaster, nMolecules-1);
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
            configuration.initializeCoordinates(box[iTemp]);
            
            setMeter(iTemp,new MeterVirial(allValueClusters[iTemp]));
//            meter[iTemp].getDataInfo().setLabel("Target/Refernce Ratio "+iTemp);
            setAccumulator(iTemp,new AccumulatorRatioAverage());

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
	public MeterVirial[] meter;
    public DataSource[] meterAccept;
    public DataSource[] meterAcceptP;
	public DataAccumulator[] accumulator;
	public DataPump[] accumulatorPump;
	public Species species;
	public ActivityIntegrate ai;
	public IntegratorMC[] integrator;
	public BoxCluster[] box;
    public ClusterAbstract[][] allValueClusters;
    public ClusterWeight[] sampleCluster;
    public MCMoveAtom[] mcMoveAtom1;
    public MCMoveBoxStep[] mcMoveRotate;
    public MCMoveBoxStep[] mcMoveMulti;
    public IntegratorPT integratorPT;

	public void setMeter(int i, MeterVirial newMeter) {
		if (accumulator[i] != null) { 
			if (accumulatorPump[i] != null) {
                integrator[i].removeIntervalAction(accumulatorPump[i]);
				accumulatorPump[i] = null;
			}
			accumulator[i] = null;
		}
		meter[i] = newMeter;
		if (meter[i] != null) {
			meter[i].setBox(box[i]);
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
		integrator[i].addIntervalAction(accumulatorPump[i]);
	}
}

