/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.models.water.PNWaterGCPM;
import etomica.models.water.SpeciesWater4PCOM;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialNonAdditiveDifference;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.virial.cluster.ClusterWeight;
import etomica.virial.cluster.ClusterWeightAbs;
import etomica.virial.wheatley.ClusterWheatleyMultibodyDerivatives;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;


/**
 * Created by Navneeth on 6/12/2017.
 */
public class ClusterWheatleyMultibodyDerivativesTest {
    int npoints = 5;
    int nder =5;
    double[] value =new double[nder+1];
    BoxCluster box;
    ClusterWheatleyMultibodyDerivatives cwmd;
    Space space = Space.getInstance(3);

    @BeforeEach
    public void setUp() {
        SpeciesGeneral speciesWater = SpeciesWater4PCOM.create(false);
        Boundary b = new BoundaryRectangularNonperiodic(space);
        final PNWaterGCPM pTarget = new PNWaterGCPM(space, b);
        MayerGeneral fTarget = new MayerGeneral(pTarget);
        PNWaterGCPM.PNWaterGCPMCached p2 = pTarget.makeCachedPairPolarization();
        PNWaterGCPM pFull = new PNWaterGCPM(space, b);
        pFull.setComponent(PNWaterGCPM.Component.INDUCTION);
        PotentialNonAdditiveDifference pnad = new PotentialNonAdditiveDifference(p2, pFull);
        MayerFunctionNonAdditiveFull fnad = new MayerFunctionNonAdditiveFull(pnad);

        cwmd = new ClusterWheatleyMultibodyDerivatives(npoints, fTarget,fnad, 1e-12, nder, true);
        ClusterWeight cl = new ClusterWeightAbs(cwmd);
        cwmd.setTemperature(1.2);

        Simulation sim = new Simulation(space);

        sim.addSpecies(speciesWater);

        box = new BoxCluster(cl, space);
        sim.addBox(box);

        box.setNMolecules(speciesWater, npoints);
        box.trialNotify();
        box.acceptNotify();
    }

    @Test
    public void testGetAllLastValues() throws Exception {

        double[] testvalue = cwmd.getAllLastValues(box);
        double[] shouldbe = new double[nder+1];
        shouldbe[0] = 0.2;

        Assertions.assertArrayEquals(shouldbe,testvalue,1e-12);

        IMoleculeList al = box.getMoleculeList();
        AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
        MoleculeChildAtomAction mcal = new MoleculeChildAtomAction(translator);
        for(int i=0; i<npoints;i++){
            translator.getTranslationVector().setX(0,i*4);
            mcal.actionPerformed(al.get(i));
        }
        box.trialNotify();
        box.acceptNotify();

        testvalue = cwmd.getAllLastValues(box);
        shouldbe = new double[] {0.1979388860563625, 0.011353771065627053, -0.06281513493500539, 0.3512268218840009, -2.0136871952313515, 12.204107895617248};

//        System.out.println(java.util.Arrays.toString(testvalue));

        Assertions.assertArrayEquals(shouldbe,testvalue,1e-12);

    }

}
