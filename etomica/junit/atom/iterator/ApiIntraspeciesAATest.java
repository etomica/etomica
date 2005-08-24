package etomica.junit.atom.iterator;

import etomica.Phase;
import etomica.Species;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.ApiIntraspeciesAA;
import etomica.junit.UnitTest;


/**
 * Unit test for ApiIntraspeciesAA
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 28, 2005 by kofke
 */
public class ApiIntraspeciesAATest extends IteratorTest {

    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 1, 6};
        int[] n2 = new int[] {1, 7, 2};
        int[] n2Tree = new int[] {3,4};
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(n0, nA0, n1, n2, n2Tree);
        AtomTreeNodeGroup rootNode = (AtomTreeNodeGroup)root.node;
        
        Species[] species = new Species[3];
        species[0] = rootNode.getDescendant(new int[] {0,0}).type.getSpecies();
        species[1] = rootNode.getDescendant(new int[] {0,1}).type.getSpecies();
        species[2] = rootNode.getDescendant(new int[] {0,2}).type.getSpecies();

        phaseTest(rootNode, species, 0);
        phaseTest(rootNode, species, 1);
        phaseTest(rootNode, species, 2);
        
        
        //test new iterator gives no iterates
        ApiIntraspeciesAA api = new ApiIntraspeciesAA(new Species[] {species[0], species[0]});
        testNoIterates(api);

        //test documented exceptions
        boolean exceptionThrown = false;
        try {
            new ApiIntraspeciesAA(new Species[] {species[0]});
        } catch(IllegalArgumentException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiIntraspeciesAA(new Species[] {species[0], species[1]});
        } catch(IllegalArgumentException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiIntraspeciesAA(new Species[] {species[0], null});
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiIntraspeciesAA((Species)null);
        } catch(NullPointerException e) {exceptionThrown = true;}
        try {
            new ApiIntraspeciesAA((Species[])null);
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);

        
    }
    
    /**
     * Performs tests on different species combinations in a particular phase.
     */
    private void phaseTest(AtomTreeNodeGroup rootNode, Species[] species, int phaseIndex) {
        speciesTestForward(rootNode, species, phaseIndex, 0);
        speciesTestForward(rootNode, species, phaseIndex, 1);
        speciesTestForward(rootNode, species, phaseIndex, 2);
    }

    /**
     * Test iteration in various directions with different targets.  Iterator constructed with
     * index of first species less than index of second.
     */
    private void speciesTestForward(AtomTreeNodeGroup rootNode, Species[] species, int phaseIndex, int species0Index) {
        ApiIntraspeciesAA api = new ApiIntraspeciesAA(species[species0Index]);
        Phase phase = rootNode.getDescendant(new int[] {phaseIndex}).node.parentPhase();
        AtomsetAction speciesTest = new SpeciesTestAction(species[species0Index], species[species0Index]);

        api.setPhase(phase);
        Atom[] molecules0 = ((AtomTreeNodeGroup)phase.getAgent(species[species0Index]).node).childList.toArray();
        
        int count = molecules0.length * (molecules0.length - 1) / 2;
        
        countTest(api, count);
        api.allAtoms(speciesTest);

        //test null phase gives no iterates
        api.setPhase(null);
        testNoIterates(api);
    }

    private class SpeciesTestAction extends AtomsetActionAdapter {
        final Species species0, species1;
        public SpeciesTestAction(Species species0, Species species1) {
            this.species0 = species0;
            this.species1 = species1;
        }
        public void actionPerformed(AtomSet atoms) {
            assertTrue(atoms.getAtom(0).type.getSpecies() == species0);
            assertTrue(atoms.getAtom(1).type.getSpecies() == species1);
        }
    }
    
}
