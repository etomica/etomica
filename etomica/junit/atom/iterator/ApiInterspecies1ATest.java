package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomManager;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.iterator.ApiInterspecies1A;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.species.Species;

/**
 * Unit test for ApiInterspecies1A
 * 
 * @author David Kofke
 *  
 */
public class ApiInterspecies1ATest extends IteratorTestAbstract {

    public void testIterator() {

        int[] n0 = new int[] { 10, 1, 0 };
        int nA0 = 5;
        int[] n1 = new int[] { 5, 1, 6 };
        int[] n2 = new int[] { 1, 7, 2 };
        int[] n2Tree = new int[] { 3, 4 };
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2,
                n2Tree);

        Species[] species = sim.getSpeciesManager().getSpecies();
        
        phaseTest(sim.getPhases()[0].getSpeciesMaster(), species);
        phaseTest(sim.getPhases()[1].getSpeciesMaster(), species);

        ApiInterspecies1A api = new ApiInterspecies1A(new Species[] {
                species[0], species[1] });

        //test new iterator gives no iterates
        testNoIterates(api);

        //test documented exceptions
        IAtom target = null;
        boolean exceptionThrown = false;
        try {
            api.setTarget(target);
        } catch (NullPointerException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            api.setTarget(null);
        } catch (NullPointerException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiInterspecies1A(new Species[] { species[0] });
        } catch (IllegalArgumentException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiInterspecies1A(new Species[] { species[0], species[0] });
        } catch (IllegalArgumentException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiInterspecies1A(new Species[] { species[0], null });
        } catch (NullPointerException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiInterspecies1A(null);
        } catch (NullPointerException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);

    }

    /**
     * Performs tests on different species combinations in a particular phase.
     */
    private void phaseTest(AtomManager atomManager, Species[] species) {
        speciesTestForward(atomManager, species[0], species[1]);
        speciesTestForward(atomManager, species[0], species[2]);
        speciesTestForward(atomManager, species[1], species[2]);
    }

    /**
     * Test iteration in various directions with different targets. Iterator
     * constructed with index of first species less than index of second.
     */
    private void speciesTestForward(AtomManager atomManager,
            Species species0, Species species1) {
        ApiInterspecies1A api = new ApiInterspecies1A(new Species[] {
                species0, species1 });
        Phase phase = atomManager.getPhase();
        AtomsetAction speciesTest = new SpeciesTestAction(
                species0, species1);
        IAtom target = null;
        IAtom targetMolecule = null;
        //test no iterates if no target
        api.setPhase(phase);
        IAtom[] molecules0 = ((AtomArrayList)phase.getAgent(species0).getChildList()).toArray();
        IAtom[] molecules1 = ((AtomArrayList)phase.getAgent(species1).getChildList()).toArray();
        int[] nMolecules = new int[] { molecules0.length, molecules1.length };
        testNoIterates(api);

        //species0 target; any direction
        target = phase.getAgent(species0).getChildList().getAtom(nMolecules[0] / 2);
        targetMolecule = target;
        api.setTarget(target);
        LinkedList list0 = testApiIterates(api, UP, targetMolecule, molecules1);
        api.allAtoms(speciesTest);

        //species0 target; up
        target = phase.getAgent(species0).getChildList().getAtom(nMolecules[0] / 2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(UP);
        testApiIterates(api, UP, targetMolecule, molecules1);
        api.allAtoms(speciesTest);

        //null direction should give previous list
        api.setDirection(null);
        LinkedList list1 = testApiIterates(api, UP, targetMolecule, molecules1);
        assertEquals(list0, list1);

        //species0 target; down
        target = phase.getAgent(species0).getChildList().getAtom(nMolecules[0] / 2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(DOWN);
        testNoIterates(api);

        //species0 leafAtom target; any direction
        if (!(species0.getMoleculeType() instanceof AtomTypeLeaf)) {
            target = ((IAtomGroup)phase.getAgent(species0).getChildList().getAtom(nMolecules[0] / 2)).getChildList().getAtom(1);
            targetMolecule = target.getParentGroup();
            api.setTarget(target);
            api.setDirection(UP);
            testApiIterates(api, UP, targetMolecule, molecules1);
            api.allAtoms(speciesTest);
        }

        //species1 target; both
        target = phase.getAgent(species1).getChildList().getAtom(nMolecules[1] / 2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(null);
        testApiIteratesSwap(api, targetMolecule, molecules0);
        api.allAtoms(speciesTest);

        //species1 target; up
        target = phase.getAgent(species1).getChildList().getAtom(nMolecules[1] / 2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(UP);
        testNoIterates(api);

        //species1 target; down
        target = phase.getAgent(species1).getChildList().getAtom(nMolecules[1] / 2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(DOWN);
        testApiIteratesSwap(api, targetMolecule, molecules0);
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
//            assertTrue(atoms.getAtom(0).type.getSpecies() == species0);
            //assertTrue(atoms.getAtom(1).type.getSpecies() == species1);
            assertTrue(atomSet.getAtom(0).getAddress() < atomSet.getAtom(1).getAddress());
        }
    }

    private final IteratorDirective.Direction UP = IteratorDirective.Direction.UP;
    private final IteratorDirective.Direction DOWN = IteratorDirective.Direction.DOWN;

}
