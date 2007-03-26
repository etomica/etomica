package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.simulation.Simulation;

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
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2,
                n2Tree);

        AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();

        //test new iterator gives no iterates
        testNoIterates(iterator);
        
        Phase[] phase = sim.getPhases();
        int[][] moleculeCount = new int[3][];
        for(int i=0; i<phase.length; i++) {
            moleculeCount[i] = new int[] {n0[i], n1[i], n2[i]};
        }

        /**
         * For each phase, check iterate atoms against full atom list, check count,
         * and for each species check count.  Check that full list is given again
         * when species is set to null
         */
        for (int i = 0; i < phase.length; i++) {
            iterator.setPhase(phase[i]);
            int count = nA0 * n0[i] + n1[i] + n2[i] * n2Tree[0] * n2Tree[1];
            LinkedList list = testIterates(iterator, phase[i].getSpeciesMaster().leafList.toArray());
            assertEquals(list.size(), phase[i].atomCount());
            assertEquals(list.size(), count);
        }

        iterator.setPhase(null);
        testNoIterates(iterator);
    }
}
