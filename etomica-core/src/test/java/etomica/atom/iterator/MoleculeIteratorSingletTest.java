/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.Molecule;
import etomica.molecule.MoleculeSetSinglet;
import etomica.molecule.iterator.MoleculeIteratorSinglet;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;

import java.util.LinkedList;


/**
 * Unit test for AtomIteratorSinglet class.
 *
 * @author David Kofke
 *
 */
public class MoleculeIteratorSingletTest extends MoleculeIteratorTestAbstract {

    private MoleculeIteratorSinglet singletIterator;
    private IMolecule testAtom1, testAtom2;
    private LinkedList list1, list2;
    
    public MoleculeIteratorSingletTest() {
        super();
    }

    public void setUp() {
        ISpecies species = new SpeciesSpheresHetero(Space3D.getInstance(), new AtomType[0]);
        singletIterator = new MoleculeIteratorSinglet();
        testAtom1 = new Molecule(species, 0);
        testAtom2 = new Molecule(species, 0);
        list1 = makeTestList(new IMoleculeList[] {new MoleculeSetSinglet(testAtom1)});
        list2 = makeTestList(new IMoleculeList[] {new MoleculeSetSinglet(testAtom2)});
    }

    public void testIterator() {
        print("starting");
        LinkedList<String> list = generalIteratorMethodTests(singletIterator);
        singletIterator.setMolecule(testAtom1);
        list = generalIteratorMethodTests(singletIterator);
        assertEquals(list,list1);
        singletIterator.setMolecule(null);
        assertNull(singletIterator.getMolecule());
        list = generalIteratorMethodTests(singletIterator);
        assertNull(singletIterator.getMolecule());
        assertTrue(list.size() == 0);
        singletIterator.setMolecule(testAtom2);
        list = generalIteratorMethodTests(singletIterator);
        assertEquals(list, list2);
        singletIterator.setMolecule(testAtom1);
        assertEquals(testAtom1, singletIterator.getMolecule());
        list = generalIteratorMethodTests(singletIterator);
        assertEquals(list, list1);
        assertEquals(testAtom1, singletIterator.getMolecule());
    }

}
