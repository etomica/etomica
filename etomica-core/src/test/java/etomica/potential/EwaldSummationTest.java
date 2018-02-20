/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.models.water.P2WaterSPCE;
import etomica.models.water.SpeciesWater3P;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Vector3D;
import etomica.units.Kelvin;
import org.junit.Before;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * This test compares the computed values against the values provided by NIST at
 * https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff.
 * NIST provides data for 4 configurations. This test checks only 1 of them. The file being compared can be changed by
 * changing the field "filenum".
 *
 * Created by Navneeth on 7/12/2017.
 */
public class EwaldSummationTest {

    private static int filenum;
    private static int[] NIST_nmol = {100,200,300,750};
    private static double[] NIST_boxl = {20,20,20,30};
    private static double[] NIST_real = {-5.58889E+05,-1.19295E+06,-1.96297E+06,-3.57226E+06};
    private static double[] NIST_fourier = {6.27009E+03, 6.03495E+03, 5.24461E+03, 7.58785E+03};
    private static double[] NIST_self = {-2.84469E+06,-5.68938E+06,-8.53407E+06,-1.42235E+07};
    private static double[] NIST_corr = {2.80999E+06, 5.61998E+06, 8.42998E+06, 1.41483E+07};

    EwaldSummation es;
    Box box;
    Simulation sim;

    public static class ChargeAgentSourceSPCE implements AtomLeafAgentManager.AgentSource<EwaldSummation.MyCharge>{
        private final Map<AtomType, EwaldSummation.MyCharge> myCharge;
        public ChargeAgentSourceSPCE(SpeciesWater3P species){
            myCharge = new HashMap<>();
            double chargeH = P2WaterSPCE.QH;
            double chargeO = P2WaterSPCE.QO;
            myCharge.put(species.getHydrogenType(), new EwaldSummation.MyCharge(chargeH));
            myCharge.put(species.getOxygenType(), new EwaldSummation.MyCharge(chargeO));
        }

        // *********************** set half(even # of particles ) as +ion, the other half -ion ***********************
        public EwaldSummation.MyCharge makeAgent(IAtom a, Box agentBox) {
            return myCharge.get(a.getType());
        }
        public void releaseAgent(EwaldSummation.MyCharge agent, IAtom atom, Box agentBox) {
            // Do nothing
        }

    }
    @Before
    public void setup() {
        Space space = Space.getInstance(3);
        sim = new Simulation(space);
        filenum = 4; // pick 1, 2, 3 or 4
        int numofmolecules = NIST_nmol[filenum - 1];
        double boxlength = NIST_boxl[filenum - 1];
        double kcut = Math.sqrt(26.999) * 2 * Math.PI / boxlength;
        double rCutRealES = 10;
        box = sim.makeBox();
        SpeciesWater3P species = new SpeciesWater3P(space, false);
        ChargeAgentSourceSPCE agentSource = new ChargeAgentSourceSPCE(species);
        AtomLeafAgentManager<EwaldSummation.MyCharge> atomAgentManager = new AtomLeafAgentManager<EwaldSummation.MyCharge>(agentSource, box);

        sim.addSpecies(species);
        box.setNMolecules(species, numofmolecules);
        box.getBoundary().setBoxSize(new Vector3D(boxlength, boxlength, boxlength));

        es = new EwaldSummation(box, atomAgentManager, space, kcut, rCutRealES);
        es.setAlpha(5.6 / boxlength);

        Configuration config = new ConfigurationResourceFile(
                String.format("etomica/potential/spce" + String.valueOf(filenum) + ".pos"),
                EwaldSummationTest.class
        );

        config.initializeCoordinates(box);


    }
    //Main method visualizes the box

    /*public static void main(String[] str) {
        EwaldSummationTest foo = new EwaldSummationTest();
        foo.setup();
        SimulationGraphic simg = new SimulationGraphic(foo.sim, foo.sim.getSpace(),null);
        simg.add(new DisplayBox(foo.sim, foo.box, foo.sim.getSpace(),null));
        simg.makeAndDisplayFrame();

    }*/

    @Test
    public void testUReal() throws Exception {
        double testval = es.uReal();
        double shouldbe = Kelvin.UNIT.toSim(NIST_real[filenum-1]);

        assertEquals(shouldbe,testval,10);

    }
    @Test
    public void testUFourier() throws Exception {
        double testval = es.uFourier();
        double shouldbe = Kelvin.UNIT.toSim(NIST_fourier[filenum-1]);

        assertEquals(shouldbe,testval,0.01);

    }
    @Test
    public void testUSelf() throws Exception {
        double testval = es.uSelf();
        double shouldbe = Kelvin.UNIT.toSim(NIST_self[filenum-1]);

        assertEquals(shouldbe,testval,100);

    }
    @Test
    public void testUBondCorr() throws Exception {
        double testval = es.uBondCorr();
        double shouldbe = Kelvin.UNIT.toSim(NIST_corr[filenum-1]);

        assertEquals(shouldbe,testval,100);

    }

}
