/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.parser;

import com.fasterxml.jackson.core.JsonProcessingException;
import etomica.action.MoleculeActionTranslateTo;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePair;
import etomica.potential.PotentialGroup;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesGeneral;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assumptions.assumeTrue;

public class ParmedStructureTest {
    private static ParmedStructure structure;

    private static final double EPSILON = 0.000001;

    @BeforeAll
    public static void setUp() throws Exception {
        assumeTrue(ParmedParser.class.getClassLoader().getResource("virtualenv/bin/parmed_json") != null);
        // this can take a while so we only do it once.
        structure = ParmedParser.parseGromacsResourceFiles("test.top", "test.gro");
    }

    @Test
    public void testGetBox() throws Exception {
        Box box = structure.getBox();
        Assertions.assertEquals(
                36935.2335047,
                box.getBoundary().volume(),
                EPSILON
        );
    }

    @Test
    public void testGetSpecies() throws Exception {
        SpeciesGeneral species = structure.getSpecies();

        Assertions.assertEquals(
                12.01078,
                species.makeMolecule().getChildList().get(0).getType().getMass(),
                EPSILON
        );

        Assertions.assertTrue(
                species.makeMolecule().getChildList().get(0).getPosition()
                .equals(new Vector3D(0.0, 0.0, 0.0))
        );
    }

    @Test
    public void testGetIntermolecularPotential() throws JsonProcessingException {
        PotentialGroup potentialGroup = structure.getIntermolecularPotential();
        SpeciesGeneral species = structure.getSpecies();
        potentialGroup.setBox(structure.getBox());

        IMolecule mol1 = species.makeMolecule();
        IMolecule mol2 = species.makeMolecule();

        MoleculeActionTranslateTo act = new MoleculeActionTranslateTo(Space3D.getInstance());
        act.setDestination(new Vector3D(5, 5, 5));
        act.actionPerformed(mol2);

        Assertions.assertEquals(
                -0.002102905446227447,
                potentialGroup.energy(new MoleculePair(
                        mol1,
                        mol2
                )),
                EPSILON
        );

    }

    @Test
    public void testGetMolecules() {
        List<IMolecule> molList = structure.getMolecules();

        Assertions.assertTrue(
                molList.get(0).getChildList().get(0).getPosition()
                .equals(new Vector3D(0, 0, 0))
        );

        Assertions.assertTrue(
                molList.get(1).getChildList().get(1).getPosition()
                .equals(new Vector3D(-0.9299999999999999, -0.56, 13.799999999999999))
        );
    }

}
