/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.AtomType;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.Molecule;
import etomica.molecule.MoleculeSetSinglet;
import etomica.molecule.iterator.MoleculeIteratorSinglet;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesSpheresHetero;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.LinkedList;

import static etomica.atom.iterator.MoleculeIteratorTestAbstract.*;
import static etomica.space3d.Space3D.getInstance;


/**
 * Unit test for AtomIteratorSinglet class.
 *
 * @author David Kofke
 */
class MoleculeIteratorSingletTest {

    private MoleculeIteratorSinglet singletIterator;
    private IMolecule testAtom1, testAtom2;
    private LinkedList list1, list2;

    public MoleculeIteratorSingletTest() {
        super();
    }

    @BeforeEach
    public void setUp() {
        ISpecies species = new SpeciesBuilder(Space3D.getInstance()).build();
        singletIterator = new MoleculeIteratorSinglet();
        testAtom1 = new Molecule(species, 0);
        testAtom2 = new Molecule(species, 0);
        list1 = makeTestList(new IMoleculeList[]{new MoleculeSetSinglet(testAtom1)});
        list2 = makeTestList(new IMoleculeList[]{new MoleculeSetSinglet(testAtom2)});
    }

    @Test
    public void testIterator() {
        print("starting");
        LinkedList<String> list = generalIteratorMethodTests(singletIterator);
        singletIterator.setMolecule(testAtom1);
        list = generalIteratorMethodTests(singletIterator);
        Assertions.assertEquals(list, list1);
        singletIterator.setMolecule(null);
        Assertions.assertNull(singletIterator.getMolecule());
        list = generalIteratorMethodTests(singletIterator);
        Assertions.assertNull(singletIterator.getMolecule());
        Assertions.assertTrue(list.size() == 0);
        singletIterator.setMolecule(testAtom2);
        list = generalIteratorMethodTests(singletIterator);
        Assertions.assertEquals(list, list2);
        singletIterator.setMolecule(testAtom1);
        Assertions.assertEquals(testAtom1, singletIterator.getMolecule());
        list = generalIteratorMethodTests(singletIterator);
        Assertions.assertEquals(list, list1);
        Assertions.assertEquals(testAtom1, singletIterator.getMolecule());
    }

}
