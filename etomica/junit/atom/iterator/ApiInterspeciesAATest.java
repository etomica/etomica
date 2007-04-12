package etomica.junit.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.iterator.ApiInterspeciesAA;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.species.Species;


/**
 * Unit test for ApiInterspeciesAA
 *
 * @author David Kofke
 *
 */
public class ApiInterspeciesAATest extends IteratorTestAbstract {

    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 1, 6};
        int[] n2 = new int[] {1, 7, 2};
        int[] n2Tree = new int[] {3,4};
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2, n2Tree);
        
        Species[] species = sim.getSpeciesManager().getSpecies();

        phaseTest(sim.getPhases()[0], species);
        phaseTest(sim.getPhases()[1], species);
        phaseTest(sim.getPhases()[2], species);
        
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
    private void phaseTest(Phase phase, Species[] species) {
        speciesTestForward(phase, species, 0, 1);
        speciesTestForward(phase, species, 0, 2);
        speciesTestForward(phase, species, 1, 2);
        speciesTestForward(phase, species, 1, 0);
        speciesTestForward(phase, species, 2, 0);
        speciesTestForward(phase, species, 2, 1);
    }

    /**
     * Test iteration in various directions with different targets.
     */
    private void speciesTestForward(Phase phase, Species[] species, int species0Index, int species1Index) {
        ApiInterspeciesAA api = new ApiInterspeciesAA(new Species[] {species[species0Index], species[species1Index]});
        AtomsetAction speciesTest = new SpeciesTestAction();

        api.setPhase(phase);
        IAtom[] molecules0 = phase.getAgent(species[species0Index]).getChildList().toArray();
        IAtom[] molecules1 = phase.getAgent(species[species1Index]).getChildList().toArray();
        
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
            assertTrue(atomSet.getAtom(0).getAddress() < atomSet.getAtom(1).getAddress());
        }
    }
    
}
