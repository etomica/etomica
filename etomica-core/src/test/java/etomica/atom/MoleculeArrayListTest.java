/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.ISpecies;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.util.Debug;
import junit.framework.TestCase;

public class MoleculeArrayListTest extends TestCase {

    protected ISpecies species = new SpeciesSpheresHetero(Space3D.getInstance(), new AtomType[0]);
    
	/*
	 * testTrimToSize()
	 */
	public void testTrimToSize() {
		final int size = 40;
		MoleculeArrayList arrayList = new MoleculeArrayList(size + 10);
		IMolecule[] listOfAtoms = new Molecule[size];
		for(int i = 0; i < size; i++) {
			listOfAtoms[i] = new Molecule(species, 0);
		    arrayList.add(listOfAtoms[i]);
		}

		arrayList.trimToSize();
		assertEquals(40, arrayList.sizeOfArray());

		for(int i = 0; i < size; i++) {
			IMolecule atom = arrayList.getMolecule(i);
			assertSame(listOfAtoms[i], atom);
		}
	}

	/*
	 * testSetGetTrimThreshold()
	 */
	public void testSetGetTrimThreshold() {
		float trimThreshold = 0.43f;
		MoleculeArrayList arrayList = new MoleculeArrayList();
		arrayList.setTrimThreshold(trimThreshold);
		float tT = arrayList.getTrimThreshold();
		assertEquals(trimThreshold, tT, .006);
	}

	/*
	 * testMaybeTrimToSize()
	 */
	public void testMaybeTrimToSize() {

		float trimThreshold = 0.5f;
		int size = 100;

		// test case at trim Threshold (should not be trimmed)
		MoleculeArrayList arrayList = new MoleculeArrayList(size);
		arrayList.setTrimThreshold(trimThreshold);
		for(int i = 0; i < size/2; i++) {
			arrayList.add(new Molecule(species, 0));
		}
        arrayList.maybeTrimToSize();
        assertEquals(size, arrayList.sizeOfArray());

        arrayList = null;

		// test case just under trim threshold (should be trimmed)
		arrayList = new MoleculeArrayList(size);
		arrayList.setTrimThreshold(trimThreshold);
		for(int i = 0; i < size/2-1; i++) {
			arrayList.add(new Molecule(species, 0));
		}
        arrayList.maybeTrimToSize();
        assertEquals(size/2 - 1, arrayList.sizeOfArray());
	}

	/*
	 * testEnsureCapacity()
	 */
	public void testEnsureCapacity() {

		int size = 20;

		// test case where capacity is already large enough.
		// Verify atoms in list are not changed.
		MoleculeArrayList arrayList = new MoleculeArrayList(size);
		IMolecule[] atomList = new IMolecule[5];
		for(int i = 0; i < 5; i++) {
			atomList[i] = new Molecule(species, 0);
			arrayList.add(atomList[i]);
		}
		arrayList.ensureCapacity(15);
		assertEquals(size, arrayList.sizeOfArray());
		for(int i = 0; i < 5; i++) {
			assertSame(atomList[i], arrayList.getMolecule(i));
		}
		arrayList = null;
		atomList = null;

		// test case where capacity needs to be increased.
		// Verify atoms in list are not changed.
		arrayList = new MoleculeArrayList(size);
		atomList = new IMolecule[5];
		for(int i = 0; i < 5; i++) {
			atomList[i] = new Molecule(species, 0);
			arrayList.add(atomList[i]);
		}
		arrayList.ensureCapacity(21); 
		assertEquals(21, arrayList.sizeOfArray());
		for(int i = 0; i < 5; i++) {
			assertSame(atomList[i], arrayList.getMolecule(i));
		}
		
	}

	/*
	 * testIsEmpty()
	 */
	public void testIsEmpty() {
		int size = 20;

		MoleculeArrayList arrayList = new MoleculeArrayList(size);	
		assertTrue(arrayList.isEmpty());

		arrayList.add(new Molecule(species, 0));
		assertFalse(arrayList.isEmpty());

		arrayList.remove(0);		
		assertTrue(arrayList.isEmpty());
	}

	/*
	 * testIndexOf()
	 */
	public void testIndexOf() {
		int size = 20;

		MoleculeArrayList arrayList = new MoleculeArrayList(size);
		IMolecule notInList = new Molecule(species, 0);
		assertEquals(-1, arrayList.indexOf(notInList));
		
		IMolecule inList = new Molecule(species, 0);
		arrayList.add(inList);
		assertEquals(0, arrayList.indexOf(inList));
		
		arrayList.remove(0);
		assertEquals(-1, arrayList.indexOf(inList));
	}

	/*
	 * testtoMoleculeArray()
	 */
	public void testtoMoleculeArray() {
		int size = 20;
		int numElems = 5;

		MoleculeArrayList arrayList = new MoleculeArrayList(size);
		
		assertEquals(0, arrayList.toMoleculeArray().length);

		IMolecule[] atomList = new IMolecule[numElems];
        for(int i = 0; i < numElems; i++) {
        	atomList[i] = new Molecule(species, 0);
        	arrayList.add(atomList[i]);
        }

        IMolecule[] aList = arrayList.toMoleculeArray();
        assertEquals(numElems, aList.length);
        for(int i = 0; i < numElems; i++) {
        	assertSame(atomList[i], aList[i]);
        }
	}

	/*
	 * testSet()
	 */
	public void testSet() {
		int size = 20;
		int numElems = 5;
        IMolecule newElem = new Molecule(species, 0);

		MoleculeArrayList arrayList = new MoleculeArrayList(size);
        IMolecule resultAtom = null;

		try {
		    resultAtom = arrayList.set(10, newElem);
		    // If the test is working properly, next line should never
		    // be executed.
		    assertNull(resultAtom);
		}
		catch (IndexOutOfBoundsException e) {
		}

		IMolecule[] atomList = new IMolecule[numElems];

        for(int i = 0; i < numElems; i++) {
        	atomList[i] = new Molecule(species, 0);
        	arrayList.add(atomList[i]);
        }

        // Replace first element in list
		try {
			resultAtom = arrayList.set(0, newElem);
			IMolecule a = arrayList.getMolecule(0);
			assertSame(a, newElem);
		}
		catch (IndexOutOfBoundsException e) {
			// Just need an assertion that will fail.
			// The generation of the exception indicates test failure.
			assertNotNull(null);
		}

        // Replace last element in list
		try {
			resultAtom = arrayList.set(4, newElem);
			IMolecule a = arrayList.getMolecule(4);
			assertSame(a, newElem);
		}
		catch (IndexOutOfBoundsException e) {
			// Just need an assertion that will fail.
			// The generation of the exception indicates test failure.
			assertNotNull(null);
		}

		try {
			resultAtom = arrayList.set(10, newElem);
			assertNull(resultAtom);
		}
		catch (IndexOutOfBoundsException e) {
		}

		try {
			resultAtom = null;
			resultAtom = arrayList.set(2, newElem);
		    assertSame(atomList[2], resultAtom);
		}
		catch (IndexOutOfBoundsException e) {
			assertNotNull(resultAtom);
		}

	}

	/*
	 * testAdd()
	 */
	public void testAdd() {
		int size = 10;
        boolean addResult;
        IMolecule[] atomList = new IMolecule[size];

		MoleculeArrayList arrayList = new MoleculeArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Molecule(species, 0);
			addResult = arrayList.add(atomList[i]);
			assertTrue(addResult);
		}

		for(int i = 0; i < size; i++) {
			assertSame(atomList[i], arrayList.getMolecule(i));
		}

		// Storage array is now full (10).  Add another atom
		IMolecule overTheTop = new Molecule(species, 0);
		addResult = arrayList.add(overTheTop);
		assertTrue(addResult);
		assertSame(overTheTop, arrayList.getMolecule(size));
		assertEquals((int)((float)size * (1.0f + MoleculeArrayList.getSizeIncreaseRatio()) + 1),
				      arrayList.sizeOfArray());

		try {
		    arrayList = new MoleculeArrayList(0);
		    addResult = arrayList.add(overTheTop);
		}
		catch(ArrayIndexOutOfBoundsException e) {
			// Just need an assertion that will fail.
			// The generation of the exception indicates test failure.
			assertNotNull(null);
		}
	}

	/*
	 * testAddAll()
	 */
	public void testAddAll() {
		int size = 10;
        IMolecule[] atomList = new IMolecule[size];
        IMolecule[] atomsetList = new IMolecule[size];

		MoleculeArrayList arrayList = new MoleculeArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Molecule(species, 0);
			arrayList.add(atomList[i]);
		}

		MoleculeArrayList atomSet = new MoleculeArrayList(size);
		for(int i = 0; i < size; i++) {
		    atomsetList[i] = new Molecule(species, 0);
			atomSet.add(atomsetList[i]);
		}
		arrayList.addAll(atomSet);
		
		for(int i = 0; i < size; i++) {
			assertSame(atomList[i], arrayList.getMolecule(i));
		}
		for(int i = 0; i < size; i++) {
			assertSame(atomsetList[i], arrayList.getMolecule(size+i));
		}
		assertEquals(23, arrayList.sizeOfArray());
	}

	/*
	 * testRemove()
	 */
	public void testRemove() {
		int size = 5;
        IMolecule[] atomList = new IMolecule[size];

		MoleculeArrayList arrayList = new MoleculeArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Molecule(species, 0);
			arrayList.add(atomList[i]);
		}

		IMolecule preRemove = null;
		try {
            preRemove = arrayList.remove(2);
		}
		catch (IndexOutOfBoundsException e) {
			System.out.println(e);
			// Exception thrown which indicates failure.
			// Fail test with any assertion that will fail.
			assertNotNull(null);
		}

        IMolecule postRemove = arrayList.getMolecule(2);

    	assertNotSame(preRemove, postRemove);

    	if (Debug.ON) {
    	    // MoleculeArrayList.getMolecule only does a range check if Debug is ON
        	try {
        		arrayList.getMolecule(size-1);
        		// If an exception is not thrown, then the test
        		// has failed.  Fail test with an assertion that
        		// will fail.
                assertNotNull(null);
    		}
        	catch (IndexOutOfBoundsException e) {
    			System.out.println(e);
    		}
    	}

	}

	/*
	 * testRemoveAndReplace()
	 */
	public void testRemoveAndReplace() {
		int size = 5;
        IMolecule[] atomList = new IMolecule[size];

		MoleculeArrayList arrayList = new MoleculeArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Molecule(species, 0);
			arrayList.add(atomList[i]);
		}

        IMolecule removeAtom = arrayList.removeAndReplace(2);
		assertSame(atomList[2], removeAtom);
		assertSame(atomList[size-1], arrayList.getMolecule(2));

        if (Debug.ON) {
            // MoleculeArrayList.getMolecule only does a range check if Debug is ON
            try {
                arrayList.getMolecule(size-1);
                // If an exception is not thrown, then the test
                // has failed.  Fail test with an assertion that
                // will fail.
                assertNotNull(null);
            }
            catch (IndexOutOfBoundsException e) {
                System.out.println(e);
            }
        }

		// Now, there are 4 items in the list.
		// from atomList, as ordered in the list :  0, 1, 4, 3 
		// So, removing the last item (index 3), should return
		// atomList[3].  Then, attempting to call getMolecule(3) should
        // NOT return an atom.
		 
		removeAtom = arrayList.removeAndReplace(3);
		assertSame(atomList[3], removeAtom);

        if (Debug.ON) {
            // MoleculeArrayList.getMolecule only does a range check if Debug is ON
    		try {
    			IMolecule atom = arrayList.getMolecule(3);
    			// Exception not thrown which indicates failure.
    			// Fail test with any assertion that will fail.
    			assertNotNull(null);
    		}
    		catch (IndexOutOfBoundsException e) {
    			System.out.println(e);
    		}
        }
	}

	/*
	 * testClear()
	 */
	public void testClear() {
		int size = 5;
        IMolecule[] atomList = new IMolecule[size];

		MoleculeArrayList arrayList = new MoleculeArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Molecule(species, 0);
			arrayList.add(atomList[i]);
		}
		assertFalse(arrayList.isEmpty());

		arrayList.clear();
		assertEquals(size, arrayList.sizeOfArray());
		assertTrue(arrayList.isEmpty());
		if (Debug.ON) {
            // MoleculeArrayList.getMolecule only does a range check if Debug is ON
    		try {
    			IMolecule atom = arrayList.getMolecule(0);
        		// If an exception is not thrown, then the test
        		// has failed.  Fail test with an assertion that
        		// will fail.
                assertNotNull(null);
    		}
        	catch (IndexOutOfBoundsException e) {
    			System.out.println(e);
    		}
		}
	}

}
