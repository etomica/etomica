package etomica.junit.atom;

import junit.framework.TestCase;
import etomica.api.IAtomLeaf;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomsetArray;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.util.Debug;

public class AtomArrayListTest extends TestCase {
	
	/*
	 * testTrimToSize()
	 */
	public void testTrimToSize() {
		ISpace space = Space.getInstance(3);
		final int size = 40;
		AtomArrayList arrayList = new AtomArrayList(size + 10);
		IAtomLeaf[] listOfAtoms = new AtomLeaf[size];
		for(int i = 0; i < size; i++) {
			listOfAtoms[i] = new AtomLeaf(space);
		    arrayList.add(listOfAtoms[i]);
		}

		arrayList.trimToSize();
		assertEquals(40, arrayList.sizeOfArray());

		for(int i = 0; i < size; i++) {
			IAtomLeaf atom = arrayList.getAtom(i);
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

		ISpace space = Space.getInstance(3);
		float trimThreshold = 0.5f;
		int size = 100;

		// test case at trim Threshold (should not be trimmed)
		AtomArrayList arrayList = new AtomArrayList(size);
		arrayList.setTrimThreshold(trimThreshold);
		for(int i = 0; i < size/2; i++) {
			arrayList.add(new AtomLeaf(space));
		}
        arrayList.maybeTrimToSize();
        assertEquals(size, arrayList.sizeOfArray());

        arrayList = null;

		// test case just under trim threshold (should be trimmed)
		arrayList = new AtomArrayList(size);
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
		AtomArrayList arrayList = new AtomArrayList(size);
		IAtomLeaf[] atomList = new IAtomLeaf[5];
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
		arrayList = new AtomArrayList(size);
		atomList = new IAtomLeaf[5];
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

		AtomArrayList arrayList = new AtomArrayList(size);	
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

		AtomArrayList arrayList = new AtomArrayList(size);
		IAtomLeaf notInList = new AtomLeaf(space);
		assertEquals(-1, arrayList.indexOf(notInList));
		
		IAtomLeaf inList = new AtomLeaf(space);
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

		AtomArrayList arrayList = new AtomArrayList(size);
		
		assertEquals(0, arrayList.toAtomLeafArray().length);

		IAtomLeaf[] atomList = new IAtomLeaf[numElems];
        for(int i = 0; i < numElems; i++) {
        	atomList[i] = new AtomLeaf(space);
        	arrayList.add(atomList[i]);
        }

        IAtomLeaf[] aList = arrayList.toAtomLeafArray();
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
        IAtomLeaf newElem = new AtomLeaf(space);

		AtomArrayList arrayList = new AtomArrayList(size);
        IAtomLeaf resultAtom = null;

		try {
		    resultAtom = arrayList.set(10, newElem);
		    // If the test is working properly, next line should never
		    // be executed.
		    assertNull(resultAtom);
		}
		catch (IndexOutOfBoundsException e) {
		}

		IAtomLeaf[] atomList = new IAtomLeaf[numElems];

        for(int i = 0; i < numElems; i++) {
        	atomList[i] = new AtomLeaf(space);
        	arrayList.add(atomList[i]);
        }

        // Replace first element in list
		try {
			resultAtom = arrayList.set(0, newElem);
			IAtomLeaf a = arrayList.getAtom(0);
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
			IAtomLeaf a = arrayList.getAtom(4);
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
		int size = 10;
        boolean addResult;
        IAtomLeaf[] atomList = new IAtomLeaf[size];

		AtomArrayList arrayList = new AtomArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new AtomLeaf(space);
			addResult = arrayList.add(atomList[i]);
			assertTrue(addResult);
		}

		for(int i = 0; i < size; i++) {
			assertSame(atomList[i], arrayList.getAtom(i));
		}

		// Storage array is now full (10).  Add another atom
		IAtomLeaf overTheTop = new AtomLeaf(space);
		addResult = arrayList.add(overTheTop);
		assertTrue(addResult);
		assertSame(overTheTop, arrayList.getAtom(size));
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
		ISpace space = Space.getInstance(3);
		int size = 10;
        IAtomLeaf[] atomList = new IAtomLeaf[size];
        IAtomLeaf[] atomsetList = new IAtomLeaf[size];

		AtomArrayList arrayList = new AtomArrayList(size);
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
		assertEquals(23, arrayList.sizeOfArray());
	}

	/*
	 * testRemove()
	 */
	public void testRemove() {
		ISpace space = Space.getInstance(3);
		int size = 5;
        IAtomLeaf[] atomList = new IAtomLeaf[size];

		AtomArrayList arrayList = new AtomArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new AtomLeaf(space);
			arrayList.add(atomList[i]);
		}

		IAtomLeaf preRemove = null;
		try {
            preRemove = arrayList.remove(2);
		}
		catch (IndexOutOfBoundsException e) {
			System.out.println(e);
			// Exception thrown which indicates failure.
			// Fail test with any assertion that will fail.
			assertNotNull(null);
		}

        IAtomLeaf postRemove = arrayList.getAtom(2);

    	assertNotSame(preRemove, postRemove);

    	if (Debug.ON) {
    	    // AtomLeafArrayList.getAtom only does a range check if Debug is ON
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

	}

	/*
	 * testRemoveAndReplace()
	 */
	public void testRemoveAndReplace() {
		ISpace space = Space.getInstance(3);
		int size = 5;
        IAtomLeaf[] atomList = new IAtomLeaf[size];

		AtomArrayList arrayList = new AtomArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new AtomLeaf(space);
			arrayList.add(atomList[i]);
		}

        IAtomLeaf removeAtom = arrayList.removeAndReplace(2);
		assertSame(atomList[2], removeAtom);
		assertSame(atomList[size-1], arrayList.getAtom(2));

        if (Debug.ON) {
            // AtomLeafArrayList.getAtom only does a range check if Debug is ON
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
    			IAtomLeaf atom = arrayList.getAtom(3);
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
		ISpace space = Space.getInstance(3);
		int size = 5;
        IAtomLeaf[] atomList = new IAtomLeaf[size];

		AtomArrayList arrayList = new AtomArrayList(size);
		for(int i = 0; i < size; i++) {
			atomList[i] = new AtomLeaf(space);
			arrayList.add(atomList[i]);
		}
		assertFalse(arrayList.isEmpty());

		arrayList.clear();
		assertEquals(size, arrayList.sizeOfArray());
		assertTrue(arrayList.isEmpty());
		if (Debug.ON) {
            // AtomLeafArrayList.getAtom only does a range check if Debug is ON
    		try {
    			IAtomLeaf atom = arrayList.getAtom(0);
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
