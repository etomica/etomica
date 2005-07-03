package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.AtomTreeNodeGroup;
import etomica.Phase;
import etomica.Species;
import etomica.SpeciesAgent;
import etomica.SpeciesRoot;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.junit.UnitTest;

/**
 * Unit test for ApiIntraspeciesAA
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jun 28, 2005 by kofke
 */
public class AtomIteratorAllMoleculesTest extends IteratorTest {

    public void testIterator() {

        int[] n0 = new int[] { 10, 1, 0, 0, 0};
        int nA0 = 5;
        int[] n1 = new int[] { 5, 1, 6, 0, 1};
        int[] n2 = new int[] { 1, 7, 2, 0, 0};
        int[] n2Tree = new int[] { 3, 4 };
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(n0, nA0, n1, n2,
                n2Tree);
        AtomTreeNodeGroup rootNode = (AtomTreeNodeGroup) root.node;

        Species[] species = new Species[3];
        species[0] = rootNode.getDescendant(new int[] { 0, 0 }).type
                .getSpecies();
        species[1] = rootNode.getDescendant(new int[] { 0, 1 }).type
                .getSpecies();
        species[2] = rootNode.getDescendant(new int[] { 0, 2 }).type
                .getSpecies();

        for(int i=0; i<n0.length; i++) {
            phaseTest(rootNode, species, i);
        }

    }

    /**
     * Performs tests on different species combinations in a particular phase.
     */
    private void phaseTest(AtomTreeNodeGroup rootNode, Species[] species, int phaseIndex) {
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules();
        Phase phase = rootNode.getDescendant(new int[] { phaseIndex }).node
                .parentPhase();

        iterator.setPhase(phase);
        
        AtomList moleculeList = new AtomList();
        for(int i=0; i<species.length; i++) {
            AtomList molecules = ((AtomTreeNodeGroup)phase.getAgent(species[i]).node).childList;
            moleculeList.addAll(new AtomIteratorListSimple(molecules));
        }
        
        LinkedList list = testIterates(iterator, moleculeList.toArray());
        assertEquals(list.size(), phase.moleculeCount());
    }
}
