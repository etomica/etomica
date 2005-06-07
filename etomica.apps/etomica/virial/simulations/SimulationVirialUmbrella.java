package etomica.virial.simulations;

import etomica.Default;
import etomica.Space;
import etomica.data.AccumulatorRatioAverage;
import etomica.potential.P2LennardJones;
import etomica.space3d.Space3D;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterWeightUmbrella;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactory;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.cluster.Standard;

/**
 * Simulation implementing the umbrella-sampling approach to evaluating a cluster
 * diagram.
 */
public class SimulationVirialUmbrella extends SimulationVirial {

	public SimulationVirialUmbrella(Space aSpace, double temperature, ClusterSum refCluster, ClusterSum[] targetClusters) {
		this(aSpace,new SpeciesFactorySpheres(),temperature,refCluster,targetClusters);
	}
	
    public SimulationVirialUmbrella(Space aSpace, SpeciesFactory speciesFactory, double temperature, ClusterSum refCluster, ClusterSum[] targetClusters) {
        super(aSpace,speciesFactory,temperature,makeUmbrellaCluster(refCluster,targetClusters),refCluster,targetClusters);
    }
    
    private static ClusterWeightUmbrella makeUmbrellaCluster(ClusterSum refSampleCluster, ClusterSum[] targetSampleClusters) {
        ClusterSum[] allSampleClusters = new ClusterSum[targetSampleClusters.length+1];
        allSampleClusters[0] = refSampleCluster;
        System.arraycopy(targetSampleClusters,0,allSampleClusters,1,targetSampleClusters.length);
        return new ClusterWeightUmbrella(allSampleClusters);
    }

	public static void main(String[] args) {
		Default.makeLJDefaults();

		final int nPoints = 5;
        double temperature = 1.3;
        double sigmaHSRef = 1.6;
        double sigmaLJ = 1.0;
        double b0 = Standard.B2HS(sigmaHSRef);
        double c0 = Standard.B3HS(sigmaHSRef);
        double d0 = Standard.B4HS(sigmaHSRef);
        Default.ATOM_SIZE = 1.0;
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+b0);
        System.out.println("B3HS: "+c0+" = "+(c0/b0/b0)+" B2HS^2");
        System.out.println("B4HS: "+d0+" = "+(d0/(b0*b0*b0))+" B2HS^3");
		
		Space3D space = new Space3D();
		
        MayerHardSphere fRef = new MayerHardSphere(space,sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(space,sigmaHSRef);
        P2LennardJones p2LJ = new P2LennardJones(space,sigmaLJ,1.0);
        System.out.println("LJ sampling");
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(space,p2LJ);
        MayerESpherical eTarget = new MayerESpherical(space,p2LJ);
        
        ClusterSum refCluster = Standard.virialCluster(nPoints, fRef, true, eRef);
        ClusterSum targetCluster = Standard.virialCluster(nPoints, fTarget, true, eTarget);
		System.out.println("B"+nPoints);

        double weightRatio = 2;
		System.out.println("Weight Ratio "+weightRatio);
		int steps = 100000000;

//		while (true) {
			SimulationVirialUmbrella sim = new SimulationVirialUmbrella(space, temperature, refCluster, 
					new ClusterSum[]{targetCluster});
			((ClusterWeightUmbrella)sim.sampleCluster).setWeightRatio(new double[] {1.0,weightRatio});
			sim.ai.setMaxSteps(steps);
//            sim.integrator.setEquilibrating(true);
			sim.ai.run();
            AccumulatorRatioAverage acc = (AccumulatorRatioAverage)sim.accumulator;
            double[][] allYourBase = (double[][])acc.getTranslator().fromArray(acc.getData());
            System.out.println("average: "+allYourBase[AccumulatorRatioAverage.RATIO.index][1]
                              +" error: "+allYourBase[AccumulatorRatioAverage.RATIO_ERROR.index][1]);
//		}
	}
}
