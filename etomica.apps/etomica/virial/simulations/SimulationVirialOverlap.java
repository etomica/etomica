package etomica.virial.simulations;

import etomica.AtomTypeGroup;
import etomica.Default;
import etomica.Integrator;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataPump;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.MCMove;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.P2LennardJones;
import etomica.space3d.Space3D;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.ConfigurationCluster;
import etomica.virial.IntegratorClusterMC;
import etomica.virial.MCMoveClusterAtom;
import etomica.virial.MCMoveClusterAtomMulti;
import etomica.virial.MCMoveClusterMolecule;
import etomica.virial.MCMoveClusterMoleculeMulti;
import etomica.virial.MCMoveClusterRotateMolecule3D;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.MeterVirial;
import etomica.virial.P0Cluster;
import etomica.virial.PhaseCluster;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

/**
 * @author kofke
 *
 * Simulation implementing the overlap-sampling approach to evaluating a cluster
 * diagram.
 */
public class SimulationVirialOverlap extends Simulation {

	public SimulationVirialOverlap(Space aSpace, SpeciesFactory speciesFactory, double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster) {
		this(aSpace,speciesFactory,temperature,new ClusterAbstract[]{refCluster,targetCluster},
                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster.makeCopy()),ClusterWeightAbs.makeWeightCluster(targetCluster.makeCopy())});
	}
	
	public SimulationVirialOverlap(Space aSpace, SpeciesFactory speciesFactory, 
            double temperature, final ClusterAbstract[] aValueClusters, final ClusterWeight[] aSampleClusters) {
		super(aSpace);

        sampleClusters = aSampleClusters;
        int nMolecules = sampleClusters[0].pointCount();
        species = speciesFactory.makeSpecies(this);
        species.setNMolecules(nMolecules);
        accumulators = new AccumulatorVirialOverlapSingleAverage[sampleClusters.length];
        accumulatorPumps = new DataPump[sampleClusters.length];
        accumulatorAAs = new IntervalActionAdapter[sampleClusters.length];
        phase = new PhaseCluster[sampleClusters.length];
        integrators = new IntegratorClusterMC[sampleClusters.length];
        meters = new MeterVirial[sampleClusters.length];
        mcMoveAtom1 = new MCMoveAtom[sampleClusters.length];
        if (nMolecules > 2) {
            mcMoveMulti = new MCMove[sampleClusters.length];
        }
        if (species.getFactory().getType() instanceof AtomTypeGroup) {
            mcMoveRotate = new MCMove[sampleClusters.length];
        }
        
        P0Cluster p0 = new P0Cluster(space);
        potentialMaster.setSpecies(p0,new Species[]{});
        
        for (int iPhase=0; iPhase<sampleClusters.length; iPhase++) {
            // integrator for iPhase samples based on iPhase cluster
            phase[iPhase] = new PhaseCluster(this,sampleClusters[iPhase]);
            
            integrators[iPhase] = new IntegratorClusterMC(potentialMaster);
            integrators[iPhase].setTemperature(temperature);
            integrators[iPhase].addPhase(phase[iPhase]);
            integrators[iPhase].setEquilibrating(false);
            if (phase[iPhase].randomMolecule().node.isLeaf()) {
                mcMoveAtom1[iPhase] = new MCMoveClusterAtom(potentialMaster);
                mcMoveAtom1[iPhase].setStepSize(1.15);
                integrators[iPhase].addMCMove(mcMoveAtom1[iPhase]);
                if (nMolecules>2) {
                    mcMoveMulti[iPhase] = new MCMoveClusterAtomMulti(potentialMaster, nMolecules-1);
                    mcMoveMulti[iPhase].setStepSize(0.41);
                    integrators[iPhase].addMCMove(mcMoveMulti[iPhase]);
                }
            }
            else {
                mcMoveAtom1[iPhase] = new MCMoveClusterMolecule(potentialMaster);
                mcMoveAtom1[iPhase].setStepSize(3.0);
                integrators[iPhase].addMCMove(mcMoveAtom1[iPhase]);
                mcMoveRotate[iPhase] = new MCMoveClusterRotateMolecule3D(potentialMaster,space);
                mcMoveRotate[iPhase].setStepSize(Math.PI);
                integrators[iPhase].addMCMove(mcMoveRotate[iPhase]);
                if (nMolecules>2) {
                    mcMoveMulti[iPhase] = new MCMoveClusterMoleculeMulti(potentialMaster, nMolecules-1);
                    mcMoveMulti[iPhase].setStepSize(0.41);
                    integrators[iPhase].addMCMove(mcMoveMulti[iPhase]);
                }
            }
            
            ConfigurationCluster configuration = new ConfigurationCluster(space);
            configuration.setPhase(phase[iPhase]);
            configuration.initializeCoordinates(phase[iPhase]);
            MeterVirial meter = new MeterVirial(new ClusterAbstract[]{aValueClusters[iPhase],aSampleClusters[1-iPhase]},integrators[iPhase]);
            setMeter(meter,iPhase);
//            meters[iPhase].getDataInfo().setLabel("Overlap/Target"+iPhase+" meter");
            AccumulatorVirialOverlapSingleAverage acc = new AccumulatorVirialOverlapSingleAverage(11);
//            acc.getDataInfo().setLabel("Overlap/Target"+iPhase+" accumulator");
            setAccumulator(acc,iPhase);
              
        }
        
        setRefPref(1,5);
        integratorOS = new IntegratorOverlap(potentialMaster, integrators, accumulators);
        integratorOS.setNumSubSteps(Default.BLOCK_SIZE);
        ai = new ActivityIntegrate(integratorOS);
        ai.setInterval(1);
        getController().addAction(ai);
		
		dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
//		dsvo.getDataInfo().setLabel("Bennet Overlap Ratio");
        integratorOS.setDSVO(dsvo);
	}

    public void setRefPref(double refPrefCenter, double span) {
        accumulators[0].setBennetParam(refPrefCenter,span);
        accumulators[1].setBennetParam(1/refPrefCenter,span);
    }
    
    public void setMeter(MeterVirial newMeter, int iPhase) {
        if (accumulators[iPhase] != null) {
            // we need a new accumulator so nuke the old one now.
            if (accumulatorPumps[iPhase] != null) {
                integrators[iPhase].removeListener(accumulatorAAs[iPhase]);
                accumulatorPumps[iPhase] = null;
            }
            accumulators[iPhase] = null;
        }
        meters[iPhase] = newMeter;
        if (newMeter != null) {
            newMeter.setPhase(phase[iPhase]);
        }
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iPhase) {
        accumulators[iPhase] = newAccumulator;
        if (accumulatorPumps[iPhase] == null) {
            accumulatorPumps[iPhase] = new DataPump(meters[iPhase],newAccumulator);
            accumulatorAAs[iPhase] = new IntervalActionAdapter(accumulatorPumps[iPhase]);
            integrators[iPhase].addListener(accumulatorAAs[iPhase]);
        }
        else {
            accumulatorPumps[iPhase].setDataSink(newAccumulator);
        }
        accumulatorAAs[iPhase].setActionInterval(1);
        if (integratorOS != null) {
            dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
//            dsvo.getDataInfo().setLabel("Bennet Overlap Ratio");
            integratorOS.setDSVO(dsvo);
        }
    }
	
    public Integrator getIntegrator() {
        return integratorOS;
    }
    
	protected DisplayPlot plot;
	public DataSourceVirialOverlap dsvo;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    protected IntervalActionAdapter[] accumulatorAAs;
    protected DataPump[] accumulatorPumps;
	protected final ClusterWeight[] sampleClusters;
    public PhaseCluster[] phase;
    protected Species species;
    public IntegratorClusterMC[] integrators;
    public MCMoveAtom[] mcMoveAtom1;
    public MCMove[] mcMoveRotate;
    public MCMove[] mcMoveMulti;
    public MeterVirial[] meters;
    public ActivityIntegrate ai;
    public IntegratorOverlap integratorOS;

	public static void main(String[] args) {
		Default.makeLJDefaults();

		int nPoints = 5;
		double temperature = 1.3; //temperature governing sampling of configurations
//		double sigmaHSRef = 1.2*sigmaLJ1B(1.0/temperature);  //diameter of reference HS system
		double sigmaHSRef = 1.6;
        double sigmaLJ = 1.0;
		double b0 = 2*Math.PI/3. * Math.pow(sigmaHSRef,3);
		Default.ATOM_SIZE = 1.0;
		System.out.println("sigmaHSRef: "+sigmaHSRef);
		System.out.println("B2HS: "+b0);
		System.out.println("B3HS: "+(5./8.*b0*b0)+" = "+(5.0/8.0)+" B2HS^2");
		System.out.println("B4HS: "+(b0*b0*b0*(219.0*Math.sqrt(2.0)/2240.0/Math.PI-89.0/280.0+4131.0/2240.0/Math.PI*Math.atan(Math.sqrt(2.0))))+" = "
				+(219.0*Math.sqrt(2.0)/2240.0/Math.PI-89.0/280.0+4131.0/2240.0/Math.PI*Math.atan(Math.sqrt(2.0)))+" B2HS^3");
		
        Space3D space = Space3D.getInstance();
		MayerHardSphere fRef = new MayerHardSphere(space,sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(space,sigmaHSRef);
        P2LennardJones p2LJ = new P2LennardJones(space,sigmaLJ,1.0);
		MayerGeneralSpherical fTarget = new MayerGeneralSpherical(space,p2LJ);
        MayerESpherical eTarget = new MayerESpherical(space,p2LJ);
		
        ClusterAbstract refCluster = Standard.virialCluster(nPoints,fRef,true,eRef,temperature);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints,fTarget,true,eTarget,temperature);

		int maxSteps = 100000;
		
        Default.BLOCK_SIZE = 1000;
//		while (true) {
			SimulationVirialOverlap sim = new SimulationVirialOverlap(space, new SpeciesFactorySpheres(), temperature, refCluster, targetCluster);
			sim.ai.setMaxSteps(maxSteps);
//            sim.integratorOS.setEquilibrating(true);
//            sim.integratorOS.setAdjustStepFreq(true);
            sim.integratorOS.setStepFreq0(0.0001);
			sim.ai.actionPerformed();
			System.out.println("average: "+sim.dsvo.getDataAsScalar()+", error: "+sim.dsvo.getError());
            DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
            System.out.println("hard sphere ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.RATIO.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.RATIO_ERROR.index)).getData()[1]);
            System.out.println("hard sphere   average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[0]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[0]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[0]);
            System.out.println("hard sphere overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);
            allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
            System.out.println("lennard jones ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.RATIO.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.RATIO_ERROR.index)).getData()[1]);
            System.out.println("lennard jones average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[0]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[0]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[0]);
            System.out.println("lennard jones overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                              +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);
//		}
		
	}//end of main
}
