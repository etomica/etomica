/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.potential.P2LennardJones;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.virial.ClusterAbstract;
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

    private static final long serialVersionUID = 1L;

    public SimulationVirialUmbrella(Space aSpace, double temperature, 
			ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
		this(aSpace,new SpeciesFactorySpheres(),temperature,refCluster,targetClusters);
	}
	
    public SimulationVirialUmbrella(Space aSpace, SpeciesFactory speciesFactory, double temperature, ClusterAbstract refCluster, ClusterAbstract[] targetClusters) {
        super(aSpace,speciesFactory,temperature,makeUmbrellaCluster(refCluster,targetClusters),refCluster,targetClusters);
    }
    
    private static ClusterWeightUmbrella makeUmbrellaCluster(ClusterAbstract refSampleCluster, ClusterAbstract[] targetSampleClusters) {
        ClusterAbstract[] allSampleClusters = new ClusterSum[targetSampleClusters.length+1];
        allSampleClusters[0] = refSampleCluster;
        System.arraycopy(targetSampleClusters,0,allSampleClusters,1,targetSampleClusters.length);
        return new ClusterWeightUmbrella(allSampleClusters);
    }

	public static void main(String[] args) {

		final int nPoints = 5;
        double temperature = 1.3;
        double sigmaHSRef = 1.6;
        double sigmaLJ = 1.0;
        double b0 = Standard.B2HS(sigmaHSRef);
        double c0 = Standard.B3HS(sigmaHSRef);
        double d0 = Standard.B4HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+b0);
        System.out.println("B3HS: "+c0+" = "+(c0/b0/b0)+" B2HS^2");
        System.out.println("B4HS: "+d0+" = "+(d0/(b0*b0*b0))+" B2HS^3");
		
		Space3D space = Space3D.getInstance();
		
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        P2LennardJones p2LJ = new P2LennardJones(space,sigmaLJ,1.0);
        System.out.println("LJ sampling");
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(p2LJ);
        MayerESpherical eTarget = new MayerESpherical(p2LJ);
        
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, true, eRef, true);
        refCluster.setTemperature(temperature);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, true, eTarget, true);
        targetCluster.setTemperature(temperature);
        
		System.out.println("B"+nPoints);

        double weightRatio = 2;
		System.out.println("Weight Ratio "+weightRatio);
		int steps = 100000000;

//		while (true) {
			SimulationVirialUmbrella sim = new SimulationVirialUmbrella(space, temperature, refCluster, 
					new ClusterAbstract[]{targetCluster});
			((ClusterWeightUmbrella)sim.sampleCluster).setWeightCoefficients(new double[] {1.0-weightRatio,weightRatio});
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
//            sim.integrator.setEquilibrating(true);
            AccumulatorRatioAverageCovariance acc = sim.accumulator;
            DataGroup allYourBase = (DataGroup)acc.getData();
            System.out.println("average: "+((DataDoubleArray)allYourBase.getData(acc.RATIO.index)).getData()[1]
                              +" error: "+((DataDoubleArray)allYourBase.getData(acc.RATIO_ERROR.index)).getData()[1]);
//		}
	}
}
