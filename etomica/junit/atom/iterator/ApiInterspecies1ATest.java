package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.ApiInterspecies1A;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.species.Species;

/**
 * Unit test for ApiInterspecies1A
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jun 28, 2005 by kofke
 */
public class ApiInterspecies1ATest extends IteratorTestAbstract {

    public void testIterator() {

        int[] n0 = new int[] { 10, 1, 0 };
        int nA0 = 5;
        int[] n1 = new int[] { 5, 1, 6 };
        int[] n2 = new int[] { 1, 7, 2 };
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
        
        phaseTest(root, species, 0);
        phaseTest(root, species, 1);

        ApiInterspecies1A api = new ApiInterspecies1A(new Species[] {
                species[0], species[1] });

        //test new iterator gives no iterates
        testNoIterates(api);

        //one species has no molecules
        api.setPhase(root.getDescendant(new int[] { 2 }).getParentPhase());
        api.setTarget(root.getDescendant(new int[] { 2, 1, 3 }));
        testNoIterates(api);
        //target not one of species
        api = new ApiInterspecies1A(new Species[] { species[1], species[2] });
        api.setPhase(root.getDescendant(new int[] { 0 }).getParentPhase());
        api.setTarget(root.getDescendant(new int[] { 0, 0, 3 }));
        testNoIterates(api);
        //target one of species but in different phase
        api = new ApiInterspecies1A(new Species[] { species[1], species[2] });
        api.setPhase(root.getDescendant(new int[] { 0 }).getParentPhase());
        api.setTarget(root.getDescendant(new int[] { 1, 1, 0 }));
        testNoIterates(api);

        //test documented exceptions
        Atom target = null;
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
    private void phaseTest(SpeciesRoot root, Species[] species,
            int phaseIndex) {
        speciesTestForward(root, species, phaseIndex, 0, 1);
        speciesTestForward(root, species, phaseIndex, 0, 2);
        speciesTestForward(root, species, phaseIndex, 1, 2);
    }

    /**
     * Test iteration in various directions with different targets. Iterator
     * constructed with index of first species less than index of second.
     */
    private void speciesTestForward(SpeciesRoot root,
            Species[] species, int phaseIndex, int species0Index,
            int species1Index) {
        ApiInterspecies1A api = new ApiInterspecies1A(new Species[] {
                species[species0Index], species[species1Index] });
        Phase phase = root.getDescendant(new int[] { phaseIndex }).getParentPhase();
        AtomsetAction speciesTest = new SpeciesTestAction(
                species[species0Index], species[species1Index]);
        Atom target = null;
        Atom targetMolecule = null;
        //test no iterates if no target
        api.setPhase(phase);
        Atom[] molecules0 = phase.getAgent(species[species0Index]).getChildList().toArray();
        Atom[] molecules1 = phase.getAgent(species[species1Index]).getChildList().toArray();
        int[] nMolecules = new int[] { molecules0.length, molecules1.length };
        testNoIterates(api);

        //species0 target; any direction
        target = root.getDescendant(new int[] { phaseIndex, species0Index,
                nMolecules[0] / 2 });
        targetMolecule = target;
        api.setTarget(target);
        LinkedList list0 = testApiIterates(api, UP, targetMolecule, molecules1);
        api.allAtoms(speciesTest);

        //species0 target; up
        target = root.getDescendant(new int[] { phaseIndex, species0Index,
                nMolecules[0] / 2 });
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
        target = root.getDescendant(new int[] { phaseIndex, species0Index,
                nMolecules[0] / 2 });
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(DOWN);
        testNoIterates(api);

        //species0 leafAtom target; any direction
        if (species0Index != 1) {
            target = root.getDescendant(new int[] { phaseIndex,
                    species0Index, nMolecules[0] / 2, 1 });
            targetMolecule = target.getParentGroup();
            api.setTarget(target);
            api.setDirection(UP);
            testApiIterates(api, UP, targetMolecule, molecules1);
            api.allAtoms(speciesTest);
        }

        //species1 target; both
        target = root.getDescendant(new int[] { phaseIndex, species1Index,
                nMolecules[1] / 2 });
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(null);
        testApiIteratesSwap(api, targetMolecule, molecules0);
        api.allAtoms(speciesTest);

        //species1 target; up
        target = root.getDescendant(new int[] { phaseIndex, species1Index,
                nMolecules[1] / 2 });
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(UP);
        testNoIterates(api);

        //species1 target; down
        target = root.getDescendant(new int[] { phaseIndex, species1Index,
                nMolecules[1] / 2 });
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(DOWN);
        testApiIteratesSwap(api, targetMolecule, molecules0);
        api.allAtoms(speciesTest);

        api.setPhase(null);
        testNoIterates(api);
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
