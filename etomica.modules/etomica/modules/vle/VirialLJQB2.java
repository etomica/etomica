package etomica.modules.vle;

import etomica.atom.iterator.ApiIntergroup;
import etomica.potential.P2LJQ;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.virial.ClusterAbstract;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactoryOrientedSpheres;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap;

/**
 * LJ simulation using Mayer sampling to evaluate cluster integrals
 */
public class VirialLJQB2 {


    public static double calcB2(double temperature, double moment) {
        if (temperature <= 0) return Double.NaN;
        int nPoints = 2;
        long steps = 1000000l;
        double sigmaHSRef = 1.5;
        moment *= moment;
        double HSB2 = Standard.B2HS(sigmaHSRef);
		
        Space space = Space3D.getInstance();
        
        PotentialGroup pTargetGroup = new PotentialGroup(2, space);
        MayerHardSphere fRef = new MayerHardSphere(space,sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(space,sigmaHSRef);
        P2LJQ pTarget = new P2LJQ(space);
        pTarget.setEpsilon(1.0);
        pTarget.setSigma(1.0);
        pTarget.setQuadrupolarMomentSquare(moment);
        pTargetGroup.addPotential(pTarget, new ApiIntergroup());
        MayerGeneral fTarget = new MayerGeneral(pTargetGroup);
        MayerEGeneral eTarget = new MayerEGeneral(pTargetGroup);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        steps /= 1000;
		
        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,new SpeciesFactoryOrientedSpheres(), temperature,refCluster,targetCluster);
        sim.integratorOS.setNumSubSteps(1000);
        // if running interactively, don't use the file
        // this will either read the refpref in from a file or run a short simulation to find it
//        sim.setRefPref(1.0082398078547523);
        sim.initRefPref(null, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(null, steps/40);
        
        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        double ratio = sim.dsvo.getDataAsScalar();
        return ratio*HSB2;
//        double error = sim.dsvo.getError();
	}
}

