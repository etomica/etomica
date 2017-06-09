/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.UnitTestUtil;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.molecule.Molecule;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.iterator.MoleculeIteratorArrayListSimple;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;

import java.util.LinkedList;

/**
 * Unit test for AtomIteratorArrayList.
 */
public class MoleculeIteratorArrayListTest extends MoleculeIteratorTestAbstract {

    public MoleculeIteratorArrayListTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }
    
    public void testListVariations() {
        ISpecies species = new SpeciesSpheresHetero(Space3D.getInstance(), new AtomType[0]);
        MoleculeIteratorArrayListSimple iterator = new MoleculeIteratorArrayListSimple();
        
        //make sure new iterator gives no iterates
        LinkedList<String> list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
        // make empty list to start
        MoleculeArrayList atomList = new MoleculeArrayList();
        iterator.setList(atomList);
        
        //add some atoms and check each time
        for(int i=0; i<=10; i++) {
            list = generalIteratorMethodTests(iterator);
            assertEquals(list.size(), i);
            atomList.add(new Molecule(species, 0));
        }
        list = generalIteratorMethodTests(iterator);
        
        //check that setList changes list
        MoleculeArrayList arrayList = new MoleculeArrayList();
        iterator.setList(arrayList);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        arrayList.add(new Molecule(species, 0));
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        
        //check handling of null list
        iterator.setList(null);
        boolean exceptionThrown = false;
        try {
            list = generalIteratorMethodTests(iterator);
        }
        catch (RuntimeException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        
    }

}

