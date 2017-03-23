package etomica.parser;

import etomica.box.Box;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresCustom;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;

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

}
