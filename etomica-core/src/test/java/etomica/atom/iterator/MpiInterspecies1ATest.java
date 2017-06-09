/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.UnitTestUtil;
import etomica.atom.MoleculesetAction;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.iterator.MpiInterspecies1A;
import etomica.potential.IteratorDirective;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;

import java.util.LinkedList;

/**
 * Unit test for ApiInterspecies1A
 * 
 * @author David Kofke
 *  
 */
public class MpiInterspecies1ATest extends MoleculeIteratorTestAbstract {

    public void testIterator() {

        int[] n0 = new int[] { 10, 1, 0 };
        int nA0 = 5;
        int[] n1 = new int[] { 5, 1, 6 };
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1);

        ISpecies[] species = new ISpecies[sim.getSpeciesCount()];
        for(int i = 0; i < sim.getSpeciesCount(); i++) {
        	species[i] = sim.getSpecies(i);
        }
        boxTest(sim.getBox(0), species);
        boxTest(sim.getBox(1), species);

        MpiInterspecies1A api = new MpiInterspecies1A(new ISpecies[] {
        		species[0], species[1]});

        //test new iterator gives no iterates
        boolean exceptionThrown = false;
        try {
            // calling reset before setting a target
            api.reset();
        } catch (NullPointerException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);

        //test documented exceptions
        IMolecule target = null;
        exceptionThrown = false;
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
            new MpiInterspecies1A(new ISpecies[] { species[0] });
        } catch (IllegalArgumentException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new MpiInterspecies1A(new ISpecies[] { species[0], species[0] });
        } catch (IllegalArgumentException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new MpiInterspecies1A(new ISpecies[] { species[0], null });
        } catch (NullPointerException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new MpiInterspecies1A(null);
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
        MpiInterspecies1A api = new MpiInterspecies1A(new ISpecies[] {
                species0, species1 });
        MoleculesetAction speciesTest = new SpeciesTestAction();
        IMolecule targetMolecule = null;
        api.setBox(box);
        IMolecule[] molecules0 = ((MoleculeArrayList)box.getMoleculeList(species0)).toMoleculeArray();
        IMolecule[] molecules1 = ((MoleculeArrayList)box.getMoleculeList(species1)).toMoleculeArray();
        int[] nMolecules = new int[] { molecules0.length, molecules1.length };

        //species0 target; any direction
        targetMolecule = box.getMoleculeList(species0).getMolecule(nMolecules[0] / 2);
        api.setTarget(targetMolecule);
        LinkedList list0 = testApiIterates(api, UP, targetMolecule, molecules1);
        MoleculeIteratorTestAbstract.allAtoms(api, speciesTest);

        //species0 target; up
        targetMolecule = box.getMoleculeList(species0).getMolecule(nMolecules[0] / 2);
        api.setTarget(targetMolecule);
        api.setDirection(UP);
        testApiIterates(api, UP, targetMolecule, molecules1);
        MoleculeIteratorTestAbstract.allAtoms(api, speciesTest);

        //null direction should give previous list
        api.setDirection(null);
        LinkedList list1 = testApiIterates(api, UP, targetMolecule, molecules1);
        assertEquals(list0, list1);

        //species0 target; down
        targetMolecule = box.getMoleculeList(species0).getMolecule(nMolecules[0] / 2);
        api.setTarget(targetMolecule);
        api.setDirection(DOWN);
        testNoIterates(api);

        //species1 target; both
        targetMolecule = box.getMoleculeList(species1).getMolecule(nMolecules[1] / 2);
        api.setTarget(targetMolecule);
        api.setDirection(null);
        testApiIteratesSwap(api, targetMolecule, molecules0);
        MoleculeIteratorTestAbstract.allAtoms(api, speciesTest);

        //species1 target; up
        targetMolecule = box.getMoleculeList(species1).getMolecule(nMolecules[1] / 2);
        api.setTarget(targetMolecule);
        api.setDirection(UP);
        testNoIterates(api);

        //species1 target; down
        targetMolecule = box.getMoleculeList(species1).getMolecule(nMolecules[1] / 2);
        api.setTarget(targetMolecule);
        api.setDirection(DOWN);
        testApiIteratesSwap(api, targetMolecule, molecules0);
        MoleculeIteratorTestAbstract.allAtoms(api, speciesTest);

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
    
    protected class SpeciesTestAction implements MoleculesetAction {

        public void actionPerformed(IMoleculeList atomSet) {
//            assertTrue(atoms.getAtom(0).type.getSpecies() == species0);
            //assertTrue(atoms.getAtom(1).type.getSpecies() == species1);
            assertTrue(atomSet.getMolecule(0).getType().getIndex() < atomSet.getMolecule(1).getType().getIndex());
        }
    }

    private final IteratorDirective.Direction UP = IteratorDirective.Direction.UP;
    private final IteratorDirective.Direction DOWN = IteratorDirective.Direction.DOWN;

}
