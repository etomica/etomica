/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.storage.DoubleStorage;
import etomica.box.storage.Tokens;
import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.models.water.P2WaterSPCE;
import etomica.models.water.SpeciesWater3P;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.RepeatedTest;
import org.junit.jupiter.api.RepetitionInfo;

import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertAll;
import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * This test compares the computed values against the values provided by NIST at
 * https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff.
 * NIST provides data for 4 configurations. This test checks only 1 of them. The file being compared can be changed by
 * changing the field "filenum".
 * <p>
 * Created by Navneeth on 7/12/2017.
 */
class EwaldSummationTest {

    private static final int[] NIST_nmol = {100, 200, 300, 750};
    private static final double[] NIST_boxl = {20, 20, 20, 30};
    private static final double[] NIST_real = {-5.58889E+05, -1.19295E+06, -1.96297E+06, -3.57226E+06};
    private static final double[] NIST_fourier = {6.27009E+03, 6.03495E+03, 5.24461E+03, 7.58785E+03};
    private static final double[] NIST_self = {-2.84469E+06, -5.68938E+06, -8.53407E+06, -1.42235E+07};
    private static final double[] NIST_corr = {2.80999E+06, 5.61998E+06, 8.42998E+06, 1.41483E+07};

    private EwaldSummation es;

    @BeforeEach
    void setup(RepetitionInfo repetitionInfo) {
        int filenum = repetitionInfo.getCurrentRepetition(); // pick 1, 2, 3 or 4
        int numofmolecules = NIST_nmol[filenum - 1];
        double boxlength = NIST_boxl[filenum - 1];
        double kcut = Math.sqrt(26.999) * 2 * Math.PI / boxlength;
        double rCutRealES = 10;
        Space space = Space.getInstance(3);
        Simulation sim = new Simulation(space);
        SpeciesGeneral species = SpeciesWater3P.create();
        sim.addSpecies(species);
        Box box = sim.makeBox();

        DoubleStorage charges = box.getAtomStorage(Tokens.doubles(new ChargeAgentSourceSPCE(species)));
        box.setNMolecules(species, numofmolecules);
        box.getBoundary().setBoxSize(new Vector3D(boxlength, boxlength, boxlength));

        es = new EwaldSummation(box, charges, space, kcut, rCutRealES);
        es.setAlpha(5.6 / boxlength);

        Configuration config = new ConfigurationResourceFile(
                String.format("spce%d.pos", filenum),
                EwaldSummationTest.class
        );

        config.initializeCoordinates(box);


    }

    @RepeatedTest(value = 4, name = "Configuration file {currentRepetition} of {totalRepetitions}")
    @DisplayName("Test Ewald Sum computed values against NIST values")
    void testEwaldSum(RepetitionInfo repetitionInfo) {
        int filenum = repetitionInfo.getCurrentRepetition();
        assertAll(
                () -> assertEquals(Kelvin.UNIT.toSim(NIST_real[filenum - 1]), es.uReal(), 10, "uReal"),
                () -> assertEquals(Kelvin.UNIT.toSim(NIST_fourier[filenum - 1]), es.uFourier(), 0.01, "uFourier"),
                () -> assertEquals(Kelvin.UNIT.toSim(NIST_self[filenum - 1]), es.uSelf(), 100, "uSelf"),
                () -> assertEquals(Kelvin.UNIT.toSim(NIST_corr[filenum - 1]), es.uBondCorr(), 100, "uBondCorr")
        );

    }

    private static class ChargeAgentSourceSPCE implements Tokens.Initializer<DoubleStorage> {
        private final Map<AtomType, Double> myCharge;

        public ChargeAgentSourceSPCE(SpeciesGeneral species) {
            myCharge = new HashMap<>();
            double chargeH = P2WaterSPCE.QH;
            double chargeO = P2WaterSPCE.QO;
            myCharge.put(species.getTypeByName("H"), chargeH);
            myCharge.put(species.getTypeByName("O"), chargeO);
        }

        // *********************** set half(even # of particles ) as +ion, the other half -ion ***********************
        @Override
        public void init(int idx, DoubleStorage storage, Box box) {
            storage.create(idx).set(myCharge.get(box.getLeafList().get(idx).getType()));
        }
    }
}
