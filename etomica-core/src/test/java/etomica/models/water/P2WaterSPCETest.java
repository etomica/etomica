package etomica.models.water;

import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.MoleculePair;
import etomica.box.Box;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by akshara on 05-10-2017.
 */
public class P2WaterSPCETest {
    private IMoleculeList molecules;
    private static final double EPSILON = 4e-7;
    @Before
    public void setUp() throws Exception {
        IMolecule mol1 = new SpeciesWater3P(Space3D.getInstance()).makeMolecule();
        IMolecule mol2 = new SpeciesWater3P(Space3D.getInstance()).makeMolecule();

        MoleculeActionTranslateTo act = new MoleculeActionTranslateTo(Space3D.getInstance());
        act.setDestination(new Vector3D(4, 4, 4));
        act.actionPerformed(mol2);

        molecules = new MoleculePair(mol1, mol2);
    }


    @Test
    public void testEnergy() throws Exception {
        P2WaterSPCE potential = new P2WaterSPCE(Space3D.getInstance());
        potential.setBox(new Box(Space3D.getInstance()));
        assertEquals(-9.315482357333167, potential.energy(molecules), EPSILON);
    }

}
