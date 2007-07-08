package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.junit.UnitTestUtil;
import etomica.box.Box;
import etomica.simulation.ISimulation;

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
        int[] n2 = new int[] { 1, 7, 2 };
        int[] n2Tree = new int[] { 3, 4 };
        ISimulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2,
                n2Tree);

        AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();

        //test new iterator gives no iterates
        testNoIterates(iterator);
        
        Box[] box = sim.getBoxs();
        int[][] moleculeCount = new int[3][];
        for(int i=0; i<box.length; i++) {
            moleculeCount[i] = new int[] {n0[i], n1[i], n2[i]};
        }

        /**
         * For each box, check iterate atoms against full atom list, check count,
         * and for each species check count.  Check that full list is given again
         * when species is set to null
         */
        for (int i = 0; i < box.length; i++) {
            iterator.setBox(box[i]);
            int count = nA0 * n0[i] + n1[i] + n2[i] * n2Tree[0] * n2Tree[1];
            LinkedList list = testIterates(iterator, ((AtomArrayList)box[i].getLeafList()).toArray());
            assertEquals(list.size(), box[i].atomCount());
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
