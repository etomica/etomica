package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.Phase;
import etomica.Species;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.junit.UnitTest;

/**
 * Unit test for AtomIteratorLeafAtoms
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on July 2, 2005 by kofke
 */
public class AtomIteratorLeafAtomsTest extends IteratorTest {

    public void testIterator() {

        int[] n0 = new int[] { 10, 1, 0 };
        int nA0 = 5;
        int[] n1 = new int[] { 5, 0, 6 };
        int[] n2 = new int[] { 1, 7, 2 };
        int[] n2Tree = new int[] { 3, 4 };
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(n0, nA0, n1, n2,
                n2Tree);
        AtomTreeNodeGroup rootNode = (AtomTreeNodeGroup) root.node;

        AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();

        //test new iterator gives no iterates
        testNoIterates(iterator);
        
        Species[] species = new Species[3];
        for(int i=0; i<species.length; i++) {
            species[i] = rootNode.getDescendant(new int[] {0,i}).type.getSpecies();
        }
        
        Phase[] phase = new Phase[3];
        int[] atomsPerMolecule = new int[] {nA0, 1, 3*4};
        int[][] moleculeCount = new int[3][];
        for(int i=0; i<phase.length; i++) {
            phase[i] = rootNode.getDescendant(new int[] { i }).node.parentPhase();
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
            LinkedList list = testIterates(iterator, phase[i].getSpeciesMaster().atomList.toArray());
            assertEquals(list.size(), phase[i].atomCount());
            assertEquals(list.size(), count);
            for(int j=0; j<species.length; j++) {
                iterator.setSpecies(species[j]);
                countTest(iterator, moleculeCount[i][j]*atomsPerMolecule[j]);
            }
            iterator.setSpecies(null);
            countTest(iterator, count);
        }

        iterator.setPhase(null);
        testNoIterates(iterator);
        iterator.setSpecies(species[0]);
        testNoIterates(iterator);
        iterator.setPhase(phase[0]);
        countTest(iterator, moleculeCount[0][0]*atomsPerMolecule[0]);
    }
}
