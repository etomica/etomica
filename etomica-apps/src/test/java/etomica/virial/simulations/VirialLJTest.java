/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate2;
import etomica.chem.elements.ElementSimple;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Virial junit test.  This just runs a simple simulation and checks that the
 * results (the Bennett alpha parameter as well as the final ratio and
 * uncertainty) are within expected limits.
 *
 * @author Andrew Schultz
 */
public class VirialLJTest {

    @Test
    public void testVirialLJ() {
        final int nPoints = 3;
        double temperature = 1;
        long steps = 1000;
        double sigmaHSRef = 1.5;

        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        Potential2Spherical pTarget = new P2LennardJones(space,1.0,1.0);
        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        MayerESpherical eTarget = new MayerESpherical(pTarget);
        ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
        targetCluster.setTemperature(temperature);
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new SpeciesSpheresMono(space, new ElementSimple("LJ")), temperature,refCluster,targetCluster);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        sim.integratorOS.setNumSubSteps(1000);
        // this will run a short simulation to find it
        sim.initRefPref(null, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(null, steps/40);
        Assertions.assertTrue(Math.abs(sim.refPref - 1.34) < 0.12, "Ref pref (alpha) within expected limits: "+sim.refPref);
        
        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integratorOS), steps);

        double[] ratioAndError = sim.dvo.getAverageAndError();
        double ratio = ratioAndError[0];
        double error = ratioAndError[1];
        System.out.println("ratio: "+ratioAndError[0]+" "+ratioAndError[1]);

        // check against expected values, 0.0604 +/- 0.0036
        Assertions.assertTrue(Math.abs(ratio - 0.0604) < 0.011, "Final ratio within expected limits: "+ratio);
        // improvements to the algorithm might lower this.  be wary of changes that raise it.
        // improvements to uncertainty estimation might alter this up or down, but it shouldn't change by much.
        Assertions.assertTrue(Math.abs(error - 0.0034) < 0.0003, "Ratio uncertainty within expected limits: "+error);
    }
}
