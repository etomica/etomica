package etomica.junit.atom;

import junit.framework.TestCase;
import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.atom.AtomArrayListTemp;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomsetArray;
import etomica.space.ISpace;
import etomica.space.Space;

public class AtomArrayListTest extends TestCase {
	
	/*
	 * testTrimToSize()
	 */
	public void testTrimToSize() {
		ISpace space = Space.getInstance(3);
		final int size = 40;
		AtomArrayListTemp arrayList = new AtomArrayListTemp(size + 10);
		IAtom[] listOfAtoms = new AtomLeaf[size];
		for(int i = 0; i < size; i++) {
			listOfAtoms[i] = new AtomLeaf(space);
		    arrayList.add(listOfAtoms[i]);
		}

		arrayList.trimToSize();
		assertEquals(40, arrayList.sizeOfArray());

		for(int i = 0; i < size; i++) {
			IAtom atom = arrayList.getAtom(i);
			assertSame(listOfAtoms[i], atom);
		}
	}

	/*
	 * testSetGetTrimThreshold()
	 */
	public void testSetGetTrimThreshold() {
		float trimThreshold = 0.43f;
		AtomArrayListTemp arrayList = new AtomArrayListTemp();
		arrayList.setTrimThreshold(trimThreshold);
		float tT = arrayList.getTrimThreshold();
		assertEquals(trimThreshold, tT, .006);
	}

	/*
	 * testMaybeTrimToSize()
	 */
	public void testMaybeTrimToSize() {

		ISpace space = Space.getInstance(3);
		float trimThreshold = 0.5f;
		int size = 100;

		// test case at trim Threshold (should not be trimmed)
		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
		arrayList.setTrimThreshold(trimThreshold);
		for(int i = 0; i < size/2; i++) {
			arrayList.add(new AtomLeaf(space));
		}
        arrayList.maybeTrimToSize();
        assertEquals(size, arrayList.sizeOfArray());

        arrayList = null;

		// test case just under trim threshold (should be trimmed)
		arrayList = new AtomArrayListTemp(size);
		arrayList.setTrimThreshold(trimThreshold);
		for(int i = 0; i < size/2-1; i++) {
			arrayList.add(new AtomLeaf(space));
		}
        arrayList.maybeTrimToSize();
        assertEquals(size/2 - 1, arrayList.sizeOfArray());
	}

	/*
	 * testEnsureCapacity()
	 */
	public void testEnsureCapacity() {

		ISpace space = Space.getInstance(3);
		int size = 20;

		// test case where capacity is already large enough.
		// Verify atoms in list are not changed.
		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
		IAtom[] atomList = new IAtom[5];
		for(int i = 0; i < 5; i++) {
			atomList[i] = new AtomLeaf(space);
			arrayList.add(atomList[i]);
		}
		arrayList.ensureCapacity(15);
		assertEquals(size, arrayList.sizeOfArray());
		for(int i = 0; i < 5; i++) {
			assertSame(atomList[i], arrayList.getAtom(i));
		}
		arrayList = null;
		atomList = null;

		// test case where capacity needs to be increased.
		// Verify atoms in list are not changed.
		arrayList = new AtomArrayListTemp(size);
		atomList = new IAtom[5];
		for(int i = 0; i < 5; i++) {
			atomList[i] = new AtomLeaf(space);
			arrayList.add(atomList[i]);
		}
		arrayList.ensureCapacity(21); 
		assertEquals(21, arrayList.sizeOfArray());
		for(int i = 0; i < 5; i++) {
			assertSame(atomList[i], arrayList.getAtom(i));
		}
		
	}

	/*
	 * testIsEmpty()
	 */
	public void testIsEmpty() {
		ISpace space = Space.getInstance(3);
		int size = 20;

		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);	
		assertTrue(arrayList.isEmpty());

		arrayList.add(new AtomLeaf(space));
		assertFalse(arrayList.isEmpty());

		arrayList.remove(0);		
		assertTrue(arrayList.isEmpty());
	}

	/*
	 * testIndexOf()
	 */
	public void testIndexOf() {
		ISpace space = Space.getInstance(3);
		int size = 20;

		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
		IAtom notInList = new AtomLeaf(space);
		assertEquals(-1, arrayList.indexOf(notInList));
		
		IAtom inList = new AtomLeaf(space);
		arrayList.add(inList);
		assertEquals(0, arrayList.indexOf(inList));
		
		arrayList.remove(0);
		assertEquals(-1, arrayList.indexOf(inList));
	}

	/*
	 * testToArray()
	 */
	public void testToArray() {
		ISpace space = Space.getInstance(3);
		int size = 20;
		int numElems = 5;

		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
		
		assertEquals(0, arrayList.toArray().length);

		IAtom[] atomList = new IAtom[numElems];
        for(int i = 0; i < numElems; i++) {
        	atomList[i] = new AtomLeaf(space);
        	arrayList.add(atomList[i]);
        }

        IAtom[] aList = arrayList.toArray();
        assertEquals(numElems, aList.length);
        for(int i = 0; i < numElems; i++) {
        	assertSame(atomList[i], aList[i]);
        }
	}

	/*
	 * testSet()
	 */
	public void testSet() {
		ISpace space = Space.getInstance(3);
		int size = 20;
		int numElems = 5;
        IAtom newElem = new AtomLeaf(space);

		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
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
        	atomList[i] = new AtomLeaf(space);
        	arrayList.add(atomList[i]);
        }

        // Replace first element in list
		try {
			resultAtom = arrayList.set(0, newElem);
			IAtom a = arrayList.getAtom(0);
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
			IAtom a = arrayList.getAtom(4);
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
		ISpace space = Space.getInstance(3);
		int size = 5;
        boolean addResult;
        IAtom[] atomList = new IAtom[size];

		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new AtomLeaf(space);
			addResult = arrayList.add(atomList[i]);
			assertTrue(addResult);
		}

		for(int i = 0; i < size; i++) {
			assertSame(atomList[i], arrayList.getAtom(i));
		}

		// Storage array is now full (5).  Add another atom
		IAtom overTheTop = new AtomLeaf(space);
		addResult = arrayList.add(overTheTop);
		assertTrue(addResult);
		assertSame(overTheTop, arrayList.getAtom(size));
		
	}

	/*
	 * testAddAll()
	 */
	public void testAddAll() {
		ISpace space = Space.getInstance(3);
		int size = 5;
        IAtom[] atomList = new IAtom[size];
        IAtom[] atomsetList = new IAtom[size];

		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new AtomLeaf(space);
			arrayList.add(atomList[i]);
		}

		AtomsetArray atomSet = new AtomsetArray(size);
		for(int i = 0; i < size; i++) {
			atomsetList[i] = new AtomLeaf(space);
		}
		atomSet.setAtoms(atomsetList);
		arrayList.addAll(atomSet);
		
		for(int i = 0; i < size; i++) {
			assertSame(atomList[i], arrayList.getAtom(i));
		}
		for(int i = 0; i < size; i++) {
			assertSame(atomsetList[i], arrayList.getAtom(size+i));
		}
	}

	/*
	 * testRemove()
	 */
	public void testRemove() {
		ISpace space = Space.getInstance(3);
		int size = 5;
        IAtom[] atomList = new IAtom[size];

		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new AtomLeaf(space);
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

        IAtom postRemove = arrayList.getAtom(2);

    	assertNotSame(preRemove, postRemove);

    	try {
    		arrayList.getAtom(size-1);
    		// If an exception is not thrown, then the test
    		// has failed.  Fail test with an assertion that
    		// will fail.
            assertNotNull(null);
		}
    	catch (IndexOutOfBoundsException e) {
			System.out.println(e);
		}

	}

	/*
	 * testRemoveAndReplace()
	 */
	public void testRemoveAndReplace() {
		ISpace space = Space.getInstance(3);
		int size = 5;
        IAtom[] atomList = new IAtom[size];

		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new AtomLeaf(space);
			arrayList.add(atomList[i]);
		}

        IAtom removeAtom = arrayList.removeAndReplace(2);
		assertSame(atomList[2], removeAtom);
		assertSame(atomList[size-1], arrayList.getAtom(2));

		try {
		    arrayList.getAtom(size-1);
			// Exception not thrown which indicates failure.
			// Fail test with any assertion that will fail.
			assertNotNull(null);
		}
		catch (IndexOutOfBoundsException e) {
			System.out.println(e);
		}

		// Now, there are 4 items in the list.
		// from atomList, as ordered in the list :  0, 1, 4, 3 
		// So, removing the last item (index 3), should return
		// atomList[3].  Then, attempting to call getAtom(3) should
        // NOT return an atom.
		 
		removeAtom = arrayList.removeAndReplace(3);
		assertSame(atomList[3], removeAtom);

		try {
			IAtom atom = arrayList.getAtom(3);
			// Exception not thrown which indicates failure.
			// Fail test with any assertion that will fail.
			assertNotNull(null);
		}
		catch (IndexOutOfBoundsException e) {
			System.out.println(e);
		}
	}

	/*
	 * testClear()
	 */
	public void testClear() {
		ISpace space = Space.getInstance(3);
		int size = 5;
        IAtom[] atomList = new IAtom[size];

		AtomArrayListTemp arrayList = new AtomArrayListTemp(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new AtomLeaf(space);
			arrayList.add(atomList[i]);
		}
		assertFalse(arrayList.isEmpty());

		arrayList.clear();
		assertEquals(size, arrayList.sizeOfArray());
		assertTrue(arrayList.isEmpty());
		try {
			IAtom atom = arrayList.getAtom(0);
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
