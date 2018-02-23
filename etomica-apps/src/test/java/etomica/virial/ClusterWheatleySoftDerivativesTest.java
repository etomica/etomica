/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.chem.elements.ElementSimple;
import etomica.potential.P2LennardJones;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * Created by Navneeth on 6/7/2017.
 */
public class ClusterWheatleySoftDerivativesTest {
    int npoints = 5;
    BoxCluster box;
    ClusterWheatleySoftDerivatives cwsd;

    @Before
    public void setup() {
        Space space = Space.getInstance(3);
        Species species = new SpeciesSpheresMono(space, new ElementSimple(""));
        P2LennardJones plj = new P2LennardJones(space);
        MayerFunction f = new MayerGeneralSpherical(plj);
        cwsd = new ClusterWheatleySoftDerivatives(npoints, f, 1e-12, 5);
        ClusterWeight cl = new ClusterWeightAbs(cwsd);
        cwsd.setTemperature(1.2);

        Simulation sim = new Simulation(space);

        box = new BoxCluster(cl, space);

        sim.addSpecies(species);
        sim.addBox(box);

        box.setNMolecules(species, npoints);
        box.trialNotify();
        box.acceptNotify();
    }
    @Test
    public void testValue() {

        double testval = cwsd.value(box);
        double shouldbe = 0.2;

        assertEquals(shouldbe,testval,1e-12);

        IAtomList al = box.getLeafList();
        for(int i=0; i<npoints;i++){
            al.get(i).getPosition().setX(0,i/2.0);
        }
        box.trialNotify();
        box.acceptNotify();

        testval = cwsd.value(box);
        shouldbe = -0.005947243342292857;

        assertEquals(shouldbe,testval,1e-12);

    }
    @Test
    public void testGetAlllastValues(){

        double[] lastval = cwsd.getAllLastValues(box);
        double[] shouldbe = new double[]{0.2,0.0,0.0,0.0,0.0,0.0};

        assertArrayEquals(lastval,shouldbe,1e-12);

        IAtomList al = box.getLeafList();
        for(int i=0; i<npoints;i++){
            al.get(i).getPosition().setX(0,i/2.0);
        }
        box.trialNotify();
        box.acceptNotify();

        lastval=cwsd.getAllLastValues(box);
        shouldbe=new double[]{-0.005947243342292857,-0.013480397126374919,-0.018697807753013534,-0.014321377449658898,-0.009692919520843135,-0.006293319547673247};

        assertArrayEquals(lastval,shouldbe,1e-12);

    }

}
