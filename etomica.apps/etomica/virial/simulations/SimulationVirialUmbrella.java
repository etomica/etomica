package etomica.virial.simulations;

import etomica.Default;
import etomica.Space;
import etomica.data.AccumulatorRatioAverage;
import etomica.potential.P2LennardJones;
import etomica.space3d.Space3D;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterWeightUmbrella;
import etomica.virial.MayerE;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.Standard;

/**
 * Simulation implementing the umbrella-sampling approach to evaluating a cluster
 * diagram.
 */
public class SimulationVirialUmbrella extends SimulationVirial {

	public SimulationVirialUmbrella(Space aSpace, double temperature, ClusterSum refCluster, ClusterSum[] targetClusters) {
		super(aSpace,temperature, makeUmbrellaCluster(refCluster,targetClusters),refCluster,targetClusters);
	}
	
    private static ClusterWeightUmbrella makeUmbrellaCluster(ClusterSum refSampleCluster, ClusterSum[] targetSampleClusters) {
        ClusterSum[] allSampleClusters = new ClusterSum[targetSampleClusters.length+1];
        allSampleClusters[0] = refSampleCluster;
        System.arraycopy(targetSampleClusters,0,allSampleClusters,1,targetSampleClusters.length);
        return new ClusterWeightUmbrella(allSampleClusters);
    }

	public static void main(String[] args) {
		Default.makeLJDefaults();
		Default.TRUNCATE_POTENTIALS = false;

		final int nPoints = 5;
        double temperature = 1.3;
        double sigmaHSRef = 1.0;
        double b0 = Standard.B2HS(sigmaHSRef);
        double c0 = Standard.C3HS(sigmaHSRef);
        double d0 = Standard.D4HS(sigmaHSRef);
        Default.ATOM_SIZE = 1.0;
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+b0);
        System.out.println("B3HS: "+c0+" = "+(c0/b0/b0)+" B2HS^2");
        System.out.println("B4HS: "+d0+" = "+(d0/(b0*b0*b0))+" B2HS^3");
		
		Space3D space = new Space3D();
		
        MayerHardSphere fRef = new MayerHardSphere(1.0);
        MayerEHardSphere eRef = new MayerEHardSphere(1.0);
        P2LennardJones p2LJ = new P2LennardJones(space,1.0,1.0,null);
        System.out.println("LJ sampling");
        MayerGeneral fTarget = new MayerGeneral(p2LJ);
        MayerE eTarget = new MayerE(p2LJ);
        
        ClusterSum refCluster = Standard.virialCluster(nPoints, fRef, true, eRef);
        ClusterSum targetCluster = Standard.virialCluster(nPoints, fTarget, true, eTarget);
		System.out.println("B"+nPoints);

        double weightRatio = 2;
		System.out.println("Weight Ratio "+weightRatio);
		int steps = 10000000;

		while (true) {
			SimulationVirialUmbrella sim = new SimulationVirialUmbrella(space, temperature, refCluster, 
					new ClusterSum[]{targetCluster});
			((ClusterWeightUmbrella)sim.sampleCluster).setWeightRatio(new double[] {1.0,weightRatio});
			sim.ai.setMaxSteps(steps);
			sim.ai.run();
            AccumulatorRatioAverage acc = (AccumulatorRatioAverage)sim.accumulator;
            double[][] allYourBase = (double[][])acc.getTranslator().fromArray(acc.getData());
            System.out.println("average: "+allYourBase[AccumulatorRatioAverage.RATIO.index][1]
                              +" error: "+allYourBase[AccumulatorRatioAverage.RATIO_ERROR.index][1]);
		}
	}
}
