/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.vle;

import etomica.action.activity.ActivityIntegrate;
import etomica.chem.elements.ElementSimple;
import etomica.potential.P2LJQ;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresRotating;
import etomica.virial.MayerGeneralAtomic;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap2;
import etomica.virial.wheatley.ClusterWheatleyHS;
import etomica.virial.wheatley.ClusterWheatleySoft;

/**
 * LJ simulation using Mayer sampling to evaluate cluster integrals
 */
public class VirialLJQB2 {


    public static double calcB2(double temperature, double moment) {
        if (temperature <= 0) return Double.NaN;
        int nPoints = 2;
        long steps = 1000000L;
        double sigmaHSRef = 1.5;
        moment *= moment;
        double HSB2 = Standard.B2HS(sigmaHSRef);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        P2LJQ pTarget = new P2LJQ(space);
        pTarget.setEpsilon(1.0);
        pTarget.setSigma(1.0);
        pTarget.setQuadrupolarMomentSquare(moment);
        MayerGeneralAtomic fTarget = new MayerGeneralAtomic(space, pTarget);
        ClusterWheatleySoft targetCluster = new ClusterWheatleySoft(nPoints, fTarget, 1e-12);
        targetCluster.setTemperature(temperature);
        ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);

        steps /= 1000;

        ISpecies species = SpeciesSpheresRotating.create(space, new ElementSimple("O"));
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species}, new int[]{nPoints}, temperature,refCluster,targetCluster);
        sim.init();
        sim.integratorOS.setNumSubSteps(1000);
        // if running interactively, don't use the file
        // this will either read the refpref in from a file or run a short simulation to find it
//        sim.setRefPref(1.0082398078547523);
        sim.initRefPref(null, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(null, steps/40);
        ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, steps);
        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(ai);

        double ratio = sim.dvo.getAverageAndError()[0];
        return ratio*HSB2;
//        double error = sim.dsvo.getError();
	}
}

