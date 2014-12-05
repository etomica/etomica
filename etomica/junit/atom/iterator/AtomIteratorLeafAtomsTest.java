/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.junit.UnitTestUtil;

/**
 * Unit test for AtomIteratorLeafAtoms
 * 
 * @author David Kofke
 *  
 */
public class AtomIteratorLeafAtomsTest extends IteratorTestAbstract {

    public void testIterator() {

        int[] n0 = new int[] { 10, 1, 0 };
        int nA0 = 5;
        int[] n1 = new int[] { 5, 0, 6 };
        ISimulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1);

        AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();

        //test new iterator gives no iterates
        testNoIterates(iterator);

        int boxCount = sim.getBoxCount();
        int[][] moleculeCount = new int[3][];
        for(int i=0; i<boxCount; i++) {
            moleculeCount[i] = new int[] {n0[i], n1[i]};
        }

        /**
         * For each box, check iterate atoms against full atom list, check count,
         * and for each species check count.  Check that full list is given again
         * when species is set to null
         */
        for (int i = 0; i < boxCount; i++) {
        	IBox box = sim.getBox(i);
            iterator.setBox(box);
            int count = nA0 * n0[i] + n1[i];
            LinkedList list = testIterates(iterator, ((AtomArrayList)box.getLeafList()).toAtomLeafArray());
            assertEquals(list.size(), box.getLeafList().getAtomCount());
            assertEquals(list.size(), count);
        }

        //test null box throws an exception
        boolean exceptionThrown = false;
        try {
            iterator.setBox(null);
        }
        catch (RuntimeException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
    }
}
