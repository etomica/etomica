package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomLeaf;
import etomica.atom.IMolecule;
import etomica.atom.iterator.ApiInterspecies1A;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.junit.UnitTestUtil;
import etomica.simulation.ISimulation;
import etomica.species.ISpecies;

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
        ISimulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1);

        ISpecies[] species = sim.getSpeciesManager().getSpecies();
        
        boxTest(sim.getBoxs()[0], species);
        boxTest(sim.getBoxs()[1], species);

        ApiInterspecies1A api = new ApiInterspecies1A(new ISpecies[] {
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
            new ApiInterspecies1A(new ISpecies[] { species[0] });
        } catch (IllegalArgumentException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiInterspecies1A(new ISpecies[] { species[0], species[0] });
        } catch (IllegalArgumentException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiInterspecies1A(new ISpecies[] { species[0], null });
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
     * Performs tests on different species combinations in a particular box.
     */
    private void boxTest(Box box, ISpecies[] species) {
        speciesTestForward(box, species[0], species[1]);
    }

    /**
     * Test iteration in various directions with different targets. Iterator
     * constructed with index of first species less than index of second.
     */
    private void speciesTestForward(Box box,
            ISpecies species0, ISpecies species1) {
        ApiInterspecies1A api = new ApiInterspecies1A(new ISpecies[] {
                species0, species1 });
        AtomsetAction speciesTest = new SpeciesTestAction();
        IAtom target = null;
        IAtom targetMolecule = null;
        //test no iterates if no target
        api.setBox(box);
        IAtom[] molecules0 = ((AtomArrayList)box.getMoleculeList(species0)).toArray();
        IAtom[] molecules1 = ((AtomArrayList)box.getMoleculeList(species1)).toArray();
        int[] nMolecules = new int[] { molecules0.length, molecules1.length };
        testNoIterates(api);

        //species0 target; any direction
        target = box.getMoleculeList(species0).getAtom(nMolecules[0] / 2);
        targetMolecule = target;
        api.setTarget(target);
        LinkedList list0 = testApiIterates(api, UP, targetMolecule, molecules1);
        api.allAtoms(speciesTest);

        //species0 target; up
        target = box.getMoleculeList(species0).getAtom(nMolecules[0] / 2);
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
        target = box.getMoleculeList(species0).getAtom(nMolecules[0] / 2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(DOWN);
        testNoIterates(api);

        //species0 leafAtom target; any direction
        target = ((IMolecule)box.getMoleculeList(species0).getAtom(nMolecules[0] / 2)).getChildList().getAtom(1);
        targetMolecule = ((IAtomLeaf)target).getParentGroup();
        api.setTarget(target);
        api.setDirection(UP);
        testApiIterates(api, UP, targetMolecule, molecules1);
        api.allAtoms(speciesTest);


        //species1 target; both
        target = box.getMoleculeList(species1).getAtom(nMolecules[1] / 2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(null);
        testApiIteratesSwap(api, targetMolecule, molecules0);
        api.allAtoms(speciesTest);

        //species1 target; up
        target = box.getMoleculeList(species1).getAtom(nMolecules[1] / 2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(UP);
        testNoIterates(api);

        //species1 target; down
        target = box.getMoleculeList(species1).getAtom(nMolecules[1] / 2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(DOWN);
        testApiIteratesSwap(api, targetMolecule, molecules0);
        api.allAtoms(speciesTest);

        //test null box throws an exception
        boolean exceptionThrown = false;
        try {
            api.setBox(null);
        }
        catch (RuntimeException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
    }

    protected class SpeciesTestAction extends AtomsetActionAdapter {

        public void actionPerformed(AtomSet atomSet) {
//            assertTrue(atoms.getAtom(0).type.getSpecies() == species0);
            //assertTrue(atoms.getAtom(1).type.getSpecies() == species1);
            assertTrue(atomSet.getAtom(0).getType().getSpecies().getIndex() < atomSet.getAtom(1).getType().getSpecies().getIndex());
        }
    }

    private final IteratorDirective.Direction UP = IteratorDirective.Direction.UP;
    private final IteratorDirective.Direction DOWN = IteratorDirective.Direction.DOWN;

}
