/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.action.MoleculeActionTranslateTo;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePair;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/**
 * Created by akshara on 05-10-2017.
 */
public class PNWaterGCPMTest {
    private IMoleculeList molecules;
    private static final double EPSILON = 7e-7;
    @BeforeEach
    public void setUp() throws Exception {
        IMolecule mol1 = SpeciesWater4PCOM.create(false).makeMolecule();
        IMolecule mol2 = SpeciesWater4PCOM.create(false).makeMolecule();

        MoleculeActionTranslateTo act = new MoleculeActionTranslateTo(Space3D.getInstance());
        act.setDestination(new Vector3D(4, 4, 4));
        act.actionPerformed(mol2);

        molecules = new MoleculePair(mol1, mol2);
    }


    @Test
    public void testEnergy() throws Exception {
        PNWaterGCPM potential = new PNWaterGCPM(Space3D.getInstance());
        potential.setBox(new Box(Space3D.getInstance()));
        Assertions.assertEquals(-14.868663567298231, potential.energy(molecules), EPSILON);
    }

}
