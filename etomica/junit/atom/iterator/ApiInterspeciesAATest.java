package etomica.junit.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.ApiInterspeciesAA;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.species.Species;


/**
 * Unit test for ApiInterspeciesAA
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 28, 2005 by kofke
 */
public class ApiInterspeciesAATest extends IteratorTestAbstract {

    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 1, 6};
        int[] n2 = new int[] {1, 7, 2};
        int[] n2Tree = new int[] {3,4};
        SpeciesRoot root = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2, n2Tree);
        AtomTreeNodeGroup rootNode = (AtomTreeNodeGroup)root.getNode();
        
        Species[] species = new Species[3];
        species[0] = rootNode.getDescendant(new int[] {0,0}).getType().getSpecies();
        species[1] = rootNode.getDescendant(new int[] {0,1}).getType().getSpecies();
        species[2] = rootNode.getDescendant(new int[] {0,2}).getType().getSpecies();

        phaseTest(rootNode, species, 0);
        phaseTest(rootNode, species, 1);
        phaseTest(rootNode, species, 2);
        
        ApiInterspeciesAA api = new ApiInterspeciesAA(new Species[] {species[0], species[1]});
        
        //test new iterator gives no iterates
        testNoIterates(api);

        //test documented exceptions
        boolean exceptionThrown = false;
        try {
            new ApiInterspeciesAA(new Species[] {species[0]});
        } catch(IllegalArgumentException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiInterspeciesAA(new Species[] {species[0], species[0]});
        } catch(IllegalArgumentException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiInterspeciesAA(new Species[] {species[0], null});
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiInterspeciesAA(null);
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);

        
    }
    
    /**
     * Performs tests on different species combinations in a particular phase.
     */
    private void phaseTest(AtomTreeNodeGroup rootNode, Species[] species, int phaseIndex) {
        speciesTestForward(rootNode, species, phaseIndex, 0, 1);
        speciesTestForward(rootNode, species, phaseIndex, 0, 2);
        speciesTestForward(rootNode, species, phaseIndex, 1, 2);
        speciesTestForward(rootNode, species, phaseIndex, 1, 0);
        speciesTestForward(rootNode, species, phaseIndex, 2, 0);
        speciesTestForward(rootNode, species, phaseIndex, 2, 1);
    }

    /**
     * Test iteration in various directions with different targets.
     */
    private void speciesTestForward(AtomTreeNodeGroup rootNode, Species[] species, int phaseIndex, int species0Index, int species1Index) {
        ApiInterspeciesAA api = new ApiInterspeciesAA(new Species[] {species[species0Index], species[species1Index]});
        Phase phase = rootNode.getDescendant(new int[] {phaseIndex}).getNode().parentPhase();
        AtomsetAction speciesTest = new SpeciesTestAction();

        api.setPhase(phase);
        Atom[] molecules0 = ((AtomTreeNodeGroup)phase.getAgent(species[species0Index]).getNode()).getChildList().toArray();
        Atom[] molecules1 = ((AtomTreeNodeGroup)phase.getAgent(species[species1Index]).getNode()).getChildList().toArray();
        
        int count = molecules0.length * molecules1.length;
        
        countTest(api, count);
        api.allAtoms(speciesTest);

        //test null phase gives no iterates
        api.setPhase(null);
        testNoIterates(api);
    }

    private class SpeciesTestAction extends AtomsetActionAdapter {
        public SpeciesTestAction() {
        }
        public void actionPerformed(AtomSet atomSet) {
            assertTrue(atomSet.getAtom(0).getNode().getAddress() < atomSet.getAtom(1).getNode().getAddress());
        }
    }
    
}
