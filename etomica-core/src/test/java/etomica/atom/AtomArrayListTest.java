/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.space3d.Space3D;
import junit.framework.TestCase;
import etomica.space.Space;
import etomica.util.Debug;

public class AtomArrayListTest extends TestCase {

	/*
	 * testTrimToSize()
	 */
	public void testTrimToSize() {
		Space space = Space.getInstance(3);
		final int size = 40;
		AtomArrayList arrayList = new AtomArrayList(size + 10);
		IAtom[] listOfAtoms = new Atom[size];
		for(int i = 0; i < size; i++) {
			listOfAtoms[i] = new Atom(space);
		    arrayList.add(listOfAtoms[i]);
		}

		arrayList.trimToSize();
		assertEquals(40, arrayList.sizeOfArray());

		for(int i = 0; i < size; i++) {
			IAtom atom = arrayList.get(i);
			assertSame(listOfAtoms[i], atom);
		}
	}

	/*
	 * testSetGetTrimThreshold()
	 */
	public void testSetGetTrimThreshold() {
		float trimThreshold = 0.43f;
		AtomArrayList arrayList = new AtomArrayList();
		arrayList.setTrimThreshold(trimThreshold);
		float tT = arrayList.getTrimThreshold();
		assertEquals(trimThreshold, tT, .006);
	}

	/*
	 * testMaybeTrimToSize()
	 */
	public void testMaybeTrimToSize() {

		Space space = Space.getInstance(3);
		float trimThreshold = 0.5f;
		int size = 100;

		// test case at trim Threshold (should not be trimmed)
		AtomArrayList arrayList = new AtomArrayList(size);
		arrayList.setTrimThreshold(trimThreshold);
		for(int i = 0; i < size/2; i++) {
			arrayList.add(new Atom(space));
		}
        arrayList.maybeTrimToSize();
        assertEquals(size, arrayList.sizeOfArray());

        arrayList = null;

		// test case just under trim threshold (should be trimmed)
		arrayList = new AtomArrayList(size);
		arrayList.setTrimThreshold(trimThreshold);
		for(int i = 0; i < size/2-1; i++) {
			arrayList.add(new Atom(space));
		}
        arrayList.maybeTrimToSize();
        assertEquals(size/2 - 1, arrayList.sizeOfArray());
	}

	/*
	 * testEnsureCapacity()
	 */
	public void testEnsureCapacity() {

		Space space = Space.getInstance(3);
		int size = 20;

		// test case where capacity is already large enough.
		// Verify atoms in list are not changed.
		AtomArrayList arrayList = new AtomArrayList(size);
		IAtom[] atomList = new IAtom[5];
		for(int i = 0; i < 5; i++) {
			atomList[i] = new Atom(space);
			arrayList.add(atomList[i]);
		}
		arrayList.ensureCapacity(15);
		assertEquals(size, arrayList.sizeOfArray());
		for(int i = 0; i < 5; i++) {
			assertSame(atomList[i], arrayList.get(i));
		}
		arrayList = null;
		atomList = null;

		// test case where capacity needs to be increased.
		// Verify atoms in list are not changed.
		arrayList = new AtomArrayList(size);
		atomList = new IAtom[5];
		for(int i = 0; i < 5; i++) {
			atomList[i] = new Atom(space);
			arrayList.add(atomList[i]);
		}
		arrayList.ensureCapacity(21);
		assertEquals(21, arrayList.sizeOfArray());
		for(int i = 0; i < 5; i++) {
			assertSame(atomList[i], arrayList.get(i));
		}

	}

	/*
	 * testIsEmpty()
	 */
	public void testIsEmpty() {
		Space space = Space.getInstance(3);
		int size = 20;

		AtomArrayList arrayList = new AtomArrayList(size);
		assertTrue(arrayList.isEmpty());

		arrayList.add(new Atom(space));
		assertFalse(arrayList.isEmpty());

		arrayList.remove(0);
		assertTrue(arrayList.isEmpty());
	}

	/*
	 * testIndexOf()
	 */
	public void testIndexOf() {
		Space space = Space.getInstance(3);
		int size = 20;

		AtomArrayList arrayList = new AtomArrayList(size);
		IAtom notInList = new Atom(space);
		assertEquals(-1, arrayList.indexOf(notInList));

		IAtom inList = new Atom(space);
		arrayList.add(inList);
		assertEquals(0, arrayList.indexOf(inList));

		arrayList.remove(0);
		assertEquals(-1, arrayList.indexOf(inList));
	}

	/*
	 * testToArray()
	 */
	public void testToArray() {
		Space space = Space.getInstance(3);
		int size = 20;
		int numElems = 5;

		AtomArrayList arrayList = new AtomArrayList(size);

		assertEquals(0, arrayList.toAtomLeafArray().length);

		IAtom[] atomList = new IAtom[numElems];
        for(int i = 0; i < numElems; i++) {
        	atomList[i] = new Atom(space);
        	arrayList.add(atomList[i]);
        }

        IAtom[] aList = arrayList.toAtomLeafArray();
        assertEquals(numElems, aList.length);
        for(int i = 0; i < numElems; i++) {
        	assertSame(atomList[i], aList[i]);
        }
	}

	/*
	 * testSet()
	 */
	public void testSet() {
		Space space = Space.getInstance(3);
		int size = 20;
		int numElems = 5;
        IAtom newElem = new Atom(space);

		AtomArrayList arrayList = new AtomArrayList(size);
        IAtom resultAtom = null;

		try {
		    resultAtom = arrayList.set(10, newElem);
		    // If the test is working properly, next line should never
		    // be executed.
		    assertNull(resultAtom);
		}
		catch (IndexOutOfBoundsException e) {
		}

		IAtom[] atomList = new IAtom[numElems];

        for(int i = 0; i < numElems; i++) {
        	atomList[i] = new Atom(space);
        	arrayList.add(atomList[i]);
        }

        // Replace first element in list
		try {
			resultAtom = arrayList.set(0, newElem);
			IAtom a = arrayList.get(0);
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
			IAtom a = arrayList.get(4);
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
		Space space = Space.getInstance(3);
		int size = 10;
        boolean addResult;
        IAtom[] atomList = new IAtom[size];

		AtomArrayList arrayList = new AtomArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Atom(space);
			addResult = arrayList.add(atomList[i]);
			assertTrue(addResult);
		}

		for(int i = 0; i < size; i++) {
			assertSame(atomList[i], arrayList.get(i));
		}

		// Storage array is now full (10).  Add another atom
		IAtom overTheTop = new Atom(space);
		addResult = arrayList.add(overTheTop);
		assertTrue(addResult);
		assertSame(overTheTop, arrayList.get(size));
		assertEquals((int)((float)size * (1.0f + AtomArrayList.getSizeIncreaseRatio()) + 1),
				      arrayList.sizeOfArray());

		try {
		    arrayList = new AtomArrayList(0);
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
		Space space = Space.getInstance(3);
		int size = 10;
        IAtom[] atomList = new IAtom[size];
        IAtom[] atomsetList = new IAtom[size];

		AtomArrayList arrayList = new AtomArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Atom(space);
			arrayList.add(atomList[i]);
		}

		AtomsetArray atomSet = new AtomsetArray(size);
		for(int i = 0; i < size; i++) {
			atomsetList[i] = new Atom(space);
		}
		atomSet.setAtoms(atomsetList);
		arrayList.addAll(atomSet);

		for(int i = 0; i < size; i++) {
			assertSame(atomList[i], arrayList.get(i));
		}
		for(int i = 0; i < size; i++) {
			assertSame(atomsetList[i], arrayList.get(size+i));
		}
		assertEquals(23, arrayList.sizeOfArray());
	}

	public void testAddAllAtomArrayList() {
		Space s = Space3D.getInstance();
		AtomArrayList l1 = new AtomArrayList();
		l1.add(new Atom(s));
		l1.add(new Atom(s));

		AtomArrayList l2 = new AtomArrayList();
        for (int i = 0; i < 20; i++) {
            l2.add(new Atom(s));
        }

        l1.addAll(l2);
        assertEquals(22, l1.size());

        for (int i = 0; i < l2.size(); i++) {
            assertSame(l2.get(i), l1.get(i + 2));
        }
	}

	/*
	 * testRemove()
	 */
	public void testRemove() {
		Space space = Space.getInstance(3);
		int size = 5;
        IAtom[] atomList = new IAtom[size];

		AtomArrayList arrayList = new AtomArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Atom(space);
			arrayList.add(atomList[i]);
		}

		IAtom preRemove = null;
		try {
            preRemove = arrayList.remove(2);
		}
		catch (IndexOutOfBoundsException e) {
			System.out.println(e);
			// Exception thrown which indicates failure.
			// Fail test with any assertion that will fail.
			assertNotNull(null);
		}

        IAtom postRemove = arrayList.get(2);

    	assertNotSame(preRemove, postRemove);

    	if (Debug.ON) {
    	    // AtomLeafArrayList.getAtom only does a range check if Debug is ON
        	try {
        		arrayList.get(size-1);
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
		Space space = Space.getInstance(3);
		int size = 5;
        IAtom[] atomList = new IAtom[size];

		AtomArrayList arrayList = new AtomArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Atom(space);
			arrayList.add(atomList[i]);
		}

        IAtom removeAtom = arrayList.removeAndReplace(2);
		assertSame(atomList[2], removeAtom);
		assertSame(atomList[size-1], arrayList.get(2));

        if (Debug.ON) {
            // AtomLeafArrayList.getAtom only does a range check if Debug is ON
            try {
                arrayList.get(size-1);
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
		// atomList[3].  Then, attempting to call getAtom(3) should
        // NOT return an atom.

		removeAtom = arrayList.removeAndReplace(3);
		assertSame(atomList[3], removeAtom);

        if (Debug.ON) {
            // AtomLeafArrayList.getAtom only does a range check if Debug is ON
    		try {
    			IAtom atom = arrayList.get(3);
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
		Space space = Space.getInstance(3);
		int size = 5;
        IAtom[] atomList = new IAtom[size];

		AtomArrayList arrayList = new AtomArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new Atom(space);
			arrayList.add(atomList[i]);
		}
		assertFalse(arrayList.isEmpty());

		arrayList.clear();
		assertEquals(size, arrayList.sizeOfArray());
		assertTrue(arrayList.isEmpty());
		if (Debug.ON) {
            // AtomLeafArrayList.getAtom only does a range check if Debug is ON
    		try {
    			IAtom atom = arrayList.get(0);
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
