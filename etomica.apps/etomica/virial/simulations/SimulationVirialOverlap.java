package etomica.virial.simulations;

import etomica.Default;
import etomica.Integrator;
import etomica.MeterAbstract;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.DataAccumulator;
import etomica.data.DataPump;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntervalActionAdapter;
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
import etomica.virial.MayerESpherical;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.MeterVirial;
import etomica.virial.P0Cluster;
import etomica.virial.PhaseCluster;
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

	public SimulationVirialOverlap(Space aSpace, double temperature, ClusterAbstract refCluster, ClusterAbstract targetCluster) {
		this(aSpace,temperature,new ClusterAbstract[]{refCluster,targetCluster},
                new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster),ClusterWeightAbs.makeWeightCluster(targetCluster)});
	}
	
	public SimulationVirialOverlap(Space aSpace, double temperature, final ClusterAbstract[] aValueClusters, final ClusterWeight[] aSampleClusters) {
		super(aSpace);

        sampleClusters = aSampleClusters;
        int nMolecules = sampleClusters[0].pointCount();
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(nMolecules);
        accumulators = new AccumulatorVirialOverlapSingleAverage[sampleClusters.length];
        accumulatorPumps = new DataPump[sampleClusters.length];
        accumulatorAAs = new IntervalActionAdapter[sampleClusters.length];
        phase = new PhaseCluster[sampleClusters.length];
        integrators = new IntegratorClusterMC[sampleClusters.length];
        meters = new MeterVirial[sampleClusters.length];
        
        P0Cluster p0 = new P0Cluster(space);
        p0.setTemperature(temperature);
        potentialMaster.setSpecies(p0,new Species[]{});
        
        for (int iPhase=0; iPhase<sampleClusters.length; iPhase++) {
            // integrator for iPhase samples based on iPhase cluster
            phase[iPhase] = new PhaseCluster(this,sampleClusters[iPhase]);
            
            integrators[iPhase] = new IntegratorClusterMC(potentialMaster);
            integrators[iPhase].setTemperature(temperature);
            integrators[iPhase].addPhase(phase[iPhase]);
            integrators[iPhase].setEquilibrating(false);
            MCMoveAtom mcMoveAtom1 = new MCMoveClusterAtom(potentialMaster);
            if (iPhase == 0) {
                mcMoveAtom1.setStepSize(0.764);
            }
            else {
                mcMoveAtom1.setStepSize(1.146);
            }
//            mcMoveAtom1.setStepSize(1.495);
            integrators[iPhase].addMCMove(mcMoveAtom1);
            if (nMolecules>2) {
                MCMoveClusterAtomMulti multiMove = new MCMoveClusterAtomMulti(potentialMaster, nMolecules-1);
                if (iPhase == 0) {
                    multiMove.setStepSize(0.261);
                }
                else {
                    multiMove.setStepSize(0.41);
                }
//                multiMove.setStepSize(0.951);
                integrators[iPhase].addMCMove(multiMove);
            }
            ConfigurationCluster configuration = new ConfigurationCluster(space);
            configuration.setPhase(phase[iPhase]);
            phase[iPhase].setConfiguration(configuration);
            MeterVirial meter = new MeterVirial(new ClusterAbstract[]{aValueClusters[iPhase],aSampleClusters[1-iPhase]},integrators[iPhase],temperature);
            setMeter(meter,iPhase);
            meters[iPhase].setLabel("Overlap/Target"+iPhase+" meter");
            AccumulatorVirialOverlapSingleAverage acc = new AccumulatorVirialOverlapSingleAverage(1);
            double refPref = 600;
            if (iPhase==1) {
                refPref = 1/refPref; // actually targetPref
            }
            acc.setBennetParam(refPref,0.1);
            acc.setLabel("Overlap/Target"+iPhase+" accumulator");
            setAccumulator(acc,iPhase);
              
        }
        
        integratorOS = new IntegratorOverlap(potentialMaster, integrators, accumulators);
        integratorOS.setNumSubSteps(Default.BLOCK_SIZE);
        integratorOS.setAdjustStepFreq(false);
        ai = new ActivityIntegrate(integratorOS);
        ai.setInterval(1);
        getController().addAction(ai);
		
		dsvo = new DataSourceVirialOverlap(accumulators[0],accumulators[1]);
		dsvo.setLabel("Bennet Overlap Ratio");
        integratorOS.setDSVO(dsvo);
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
            newMeter.setPhase(new Phase[]{phase[iPhase]});
        }
    }

    public void setAccumulator(AccumulatorVirialOverlapSingleAverage newAccumulator, int iPhase) {
        accumulators[iPhase] = newAccumulator;
        if (accumulatorPumps[iPhase] == null) {
            accumulatorPumps[iPhase] = new DataPump(meters[iPhase],new DataAccumulator[] {newAccumulator});
            accumulatorAAs[iPhase] = new IntervalActionAdapter(accumulatorPumps[iPhase]);
        }
        else {
            accumulatorPumps[iPhase].setDataSinks(new AccumulatorRatioAverage[] {newAccumulator});
        }
        integrators[iPhase].addListener(accumulatorAAs[iPhase]);
        accumulatorAAs[iPhase].setActionInterval(1);
    }
	
    public Integrator getIntegrator() {
        return integratorOS;
    }
    
	protected DisplayPlot plot;
	protected DataSourceVirialOverlap dsvo;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    protected IntervalActionAdapter[] accumulatorAAs;
    protected DataPump[] accumulatorPumps;
	protected final ClusterWeight[] sampleClusters;
    protected PhaseCluster[] phase;
    protected Species species;
    protected IntegratorClusterMC[] integrators;
    protected MeterAbstract[] meters;
    protected ActivityIntegrate ai;
    protected IntegratorOverlap integratorOS;

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
		
        Space3D space = new Space3D();
		MayerHardSphere fRef = new MayerHardSphere(space,sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(space,sigmaHSRef);
        P2LennardJones p2LJ = new P2LennardJones(space,sigmaLJ,1.0);
		MayerGeneralSpherical fTarget = new MayerGeneralSpherical(space,p2LJ);
        MayerESpherical eTarget = new MayerESpherical(space,p2LJ);
		
        ClusterAbstract refCluster = Standard.virialCluster(nPoints,fRef,true,eRef);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints,fTarget,true,eTarget);

		int maxSteps = 100000;
		
        Default.BLOCK_SIZE = 1000;
//		while (true) {
			SimulationVirialOverlap sim = new SimulationVirialOverlap(space, temperature, refCluster, targetCluster);
			sim.ai.setMaxSteps(maxSteps);
//            sim.integratorOS.setEquilibrating(true);
//            sim.integratorOS.setAdjustStepFreq(true);
            sim.integratorOS.setStepFreq0(0.0001);
			sim.ai.actionPerformed();
			System.out.println("average: "+sim.dsvo.getData()[0]+", error: "+sim.dsvo.getError());
            double[][] allYourBase = (double[][])sim.accumulators[0].getTranslator().fromArray(sim.accumulators[0].getData(sim.dsvo.minDiffLocation()));
            System.out.println("hard sphere ratio average: "+allYourBase[AccumulatorRatioAverage.RATIO.index][1]
                              +" error: "+allYourBase[AccumulatorRatioAverage.RATIO_ERROR.index][1]);
            System.out.println("hard sphere   average: "+allYourBase[AccumulatorAverage.AVERAGE.index][0]
                              +" stdev: "+allYourBase[AccumulatorAverage.STANDARD_DEVIATION.index][0]
                              +" error: "+allYourBase[AccumulatorAverage.ERROR.index][0]);
            System.out.println("hard sphere overlap average: "+allYourBase[AccumulatorAverage.AVERAGE.index][1]
                              +" stdev: "+allYourBase[AccumulatorAverage.STANDARD_DEVIATION.index][1]
                              +" error: "+allYourBase[AccumulatorAverage.ERROR.index][1]);
            allYourBase = (double[][])sim.accumulators[1].getTranslator().fromArray(sim.accumulators[1].getData(sim.dsvo.minDiffLocation()));
            System.out.println("lennard jones ratio average: "+allYourBase[AccumulatorRatioAverage.RATIO.index][1]
                              +" error: "+allYourBase[AccumulatorRatioAverage.RATIO_ERROR.index][1]);
            System.out.println("lennard jones average: "+allYourBase[AccumulatorAverage.AVERAGE.index][0]
                              +" stdev: "+allYourBase[AccumulatorAverage.STANDARD_DEVIATION.index][0]
                              +" error: "+allYourBase[AccumulatorAverage.ERROR.index][0]);
            System.out.println("lennard jones overlap average: "+allYourBase[AccumulatorAverage.AVERAGE.index][1]
                              +" stdev: "+allYourBase[AccumulatorAverage.STANDARD_DEVIATION.index][1]
                              +" error: "+allYourBase[AccumulatorAverage.ERROR.index][1]);
//		}
		
	}//end of main
}
