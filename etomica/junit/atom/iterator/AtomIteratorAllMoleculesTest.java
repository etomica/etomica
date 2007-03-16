package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.SpeciesMaster;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.species.Species;

/**
 * Unit test for ApiIntraspeciesAA
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jun 28, 2005 by kofke
 */
public class AtomIteratorAllMoleculesTest extends IteratorTestAbstract {

    public void testIterator() {

        int[] n0 = new int[] { 10, 1, 0, 0, 0};
        int nA0 = 5;
        int[] n1 = new int[] { 5, 1, 6, 0, 1};
        int[] n2 = new int[] { 1, 7, 2, 0, 0};
        int[] n2Tree = new int[] { 3, 4 };
        SpeciesRoot root = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2,
                n2Tree);

        Species[] species = new Species[3];
        species[0] = root.getDescendant(new int[] { 0, 0 }).getType()
                .getSpecies();
        species[1] = root.getDescendant(new int[] { 0, 1 }).getType()
                .getSpecies();
        species[2] = root.getDescendant(new int[] { 0, 2 }).getType()
                .getSpecies();

        for(int i=0; i<n0.length; i++) {
            phaseTest(root, species, i);
        }

    }

    /**
     * Performs tests on different species combinations in a particular phase.
     */
    private void phaseTest(SpeciesRoot root, Species[] species, int phaseIndex) {
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules();
        Phase phase = ((SpeciesMaster)root.getChildList().get(phaseIndex)).getPhase();

        iterator.setPhase(phase);
        
        AtomArrayList moleculeList = new AtomArrayList();
        for(int i=0; i<species.length; i++) {
            AtomArrayList molecules = phase.getAgent(species[i]).getChildList();
            for (int j=0; j<molecules.size(); j++) {
                moleculeList.add(molecules.get(j));
            }
        }
        
        LinkedList list = testIterates(iterator, moleculeList.toArray());
        assertEquals(list.size(), phase.moleculeCount());
    }
}
