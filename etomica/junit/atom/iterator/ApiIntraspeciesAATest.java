package etomica.junit.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.iterator.ApiIntraspeciesAA;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.species.Species;


/**
 * Unit test for ApiIntraspeciesAA
 *
 * @author David Kofke
 *
 */

public class ApiIntraspeciesAATest extends IteratorTestAbstract {

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
    private void phaseTest(Phase phase, Species[] species) {
        speciesTestForward(phase, species, 0);
        speciesTestForward(phase, species, 1);
        speciesTestForward(phase, species, 2);
    }

    /**
     * Test iteration in various directions with different targets.  Iterator constructed with
     * index of first species less than index of second.
     */
    private void speciesTestForward(Phase phase, Species[] species, int species0Index) {
        ApiIntraspeciesAA api = new ApiIntraspeciesAA(species[species0Index]);
        AtomsetAction speciesTest = new SpeciesTestAction(species[species0Index], species[species0Index]);

        api.setPhase(phase);
        IAtom[] molecules0 = ((AtomArrayList)phase.getAgent(species[species0Index]).getChildList()).toArray();
        
        int count = molecules0.length * (molecules0.length - 1) / 2;
        
        countTest(api, count);
        api.allAtoms(speciesTest);

        //test null phase throws an exception
        boolean exceptionThrown = false;
        try {
            api.setPhase(null);
        }
        catch (RuntimeException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
    }

    private class SpeciesTestAction extends AtomsetActionAdapter {
        final Species species0, species1;
        public SpeciesTestAction(Species species0, Species species1) {
            this.species0 = species0;
            this.species1 = species1;
        }
        public void actionPerformed(AtomSet atomSet) {
            assertTrue(atomSet.getAtom(0).getType().getSpecies() == species0);
            assertTrue(atomSet.getAtom(1).getType().getSpecies() == species1);
        }
    }
    
}
