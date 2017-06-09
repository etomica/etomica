package etomica.parser;

import com.fasterxml.jackson.core.JsonProcessingException;
import etomica.action.MoleculeActionTranslateTo;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePair;
import etomica.potential.PotentialGroup;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresCustom;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class ParmedStructureTest {
    private static ParmedStructure structure;

    private static final double EPSILON = 0.000001;

    @BeforeClass
    public static void setUp() throws Exception {
        // this can take a while so we only do it once.
        structure = ParmedParser.parseGromacsResourceFiles("test.top", "test.gro");
    }

    @Test
    public void testGetBox() throws Exception {
        Box box = structure.getBox();
        assertEquals(
                36935.2335047,
                box.getBoundary().volume(),
                EPSILON
        );
    }

    @Test
    public void testGetSpecies() throws Exception {
        SpeciesSpheresCustom species = structure.getSpecies();

        assertEquals(
                12.01078,
                species.makeMolecule().getChildList().getAtom(0).getType().getMass(),
                EPSILON
        );

        assertTrue(
                species.makeMolecule().getChildList().getAtom(0).getPosition()
                .equals(new Vector3D(0.0, 0.0, 0.0))
        );
    }

    @Test
    public void testGetIntermolecularPotential() throws JsonProcessingException {
        PotentialGroup potentialGroup = structure.getIntermolecularPotential();
        SpeciesSpheresCustom species = structure.getSpecies();
        potentialGroup.setBox(structure.getBox());

        IMolecule mol1 = species.makeMolecule();
        IMolecule mol2 = species.makeMolecule();

        MoleculeActionTranslateTo act = new MoleculeActionTranslateTo(Space3D.getInstance());
        act.setDestination(new Vector3D(5, 5, 5));
        act.actionPerformed(mol2);

        assertEquals(
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

        assertTrue(
                molList.get(0).getChildList().getAtom(0).getPosition()
                .equals(new Vector3D(0, 0, 0))
        );

        assertTrue(
                molList.get(1).getChildList().getAtom(1).getPosition()
                .equals(new Vector3D(-0.9299999999999999, -0.56, 13.799999999999999))
        );
    }

}
