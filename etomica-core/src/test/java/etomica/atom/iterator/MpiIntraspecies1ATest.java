/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.MoleculesetAction;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.iterator.MpiIntraspecies1A;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;
import org.junit.jupiter.api.Test;

import static etomica.UnitTestUtil.makeStandardSpeciesTree;
import static etomica.atom.iterator.MoleculeIteratorTestAbstract.*;
import static etomica.potential.IteratorDirective.Direction;
import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * Unit test for Apiintraspecies1A
 *
 * @author David Kofke
 */
class MpiIntraspecies1ATest {

    @Test
    public void testIterator() {

        int[] n0 = new int[]{10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[]{5, 1, 6};
        Simulation sim = makeStandardSpeciesTree(n0, nA0, n1);

        ISpecies[] species = new ISpecies[sim.getSpeciesCount()];
        for (int i = 0; i < sim.getSpeciesCount(); i++) {
            species[i] = sim.getSpecies(i);
        }
        boxTest(sim.getBox(0), species);
        boxTest(sim.getBox(1), species);

        MpiIntraspecies1A api = new MpiIntraspecies1A(species[0]);

        //test documented exceptions
        boolean exceptionThrown = false;
        try {
            // need to call setTarget before calling reset
            api.reset();
        } catch (RuntimeException e) {
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
            api.setBox(null);
        } catch (RuntimeException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new MpiIntraspecies1A(null);
        } catch (NullPointerException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);


    }

    /**
     * Performs tests on different species combinations in a particular box.
     */
    private void boxTest(Box box, ISpecies[] species) {
        speciesTestForward(box, species[0]);
        speciesTestForward(box, species[1]);
    }

    /**
     * Test iteration in various directions with different targets.  Iterator constructed with index of first species
     * less than index of second.
     */
    private void speciesTestForward(Box box, ISpecies species) {

        MpiIntraspecies1A api = new MpiIntraspecies1A(species);
        MoleculesetAction speciesTest = new SpeciesTestAction(species, species);
        IMolecule target = null;
        //test no iterates if no target
        api.setBox(box);

        IMolecule[] molecules0 = ((MoleculeArrayList) box.getMoleculeList(species)).toMoleculeArray();
        int[] nMolecules = new int[]{molecules0.length};

        if (nMolecules[0] == 0) return;

        //species0 target; any direction

        target = box.getMoleculeList(species).get(nMolecules[0] / 2);
        api.setTarget(target);
        testApiIterates(api, target, upMolecules(target, molecules0), dnMolecules(target, molecules0));
        allAtoms(api, speciesTest);

        //species0 target; up
        target = box.getMoleculeList(species).get(nMolecules[0] / 2);
        api.setTarget(target);
        api.setDirection(UP);
        testApiIterates(api, UP, target, upMolecules(target, molecules0));
        allAtoms(api, speciesTest);

        //species0 target; down
        target = box.getMoleculeList(species).get(nMolecules[0] / 2);
        api.setTarget(target);
        api.setDirection(DOWN);
        testApiIteratesSwap(api, target, dnMolecules(target, molecules0));
    }

    private class SpeciesTestAction implements MoleculesetAction {
        final ISpecies species0, species1;

        public SpeciesTestAction(ISpecies species0, ISpecies species1) {
            this.species0 = species0;
            this.species1 = species1;
        }

        public void actionPerformed(IMoleculeList atomSet) {
            assertTrue(atomSet.get(0).getType() == species0);
            assertTrue(atomSet.get(1).getType() == species1);
        }
    }

    private IMolecule[] upMolecules(IMolecule target, IMolecule[] list) {
        int i;
        for (i = 0; i < list.length; i++) {
            if (list[i] == target) break;
        }
        IMolecule[] atoms = new IMolecule[list.length - i - 1];
        for (int j = 0; j < atoms.length; j++) {
            atoms[j] = list[i + j + 1];
        }
        return atoms;
    }

    private IMolecule[] dnMolecules(IMolecule target, IMolecule[] list) {
        int i;
        for (i = 0; i < list.length; i++) {
            if (list[i] == target) break;
        }
        IMolecule[] atoms = new IMolecule[i];
        for (int j = 0; j < i; j++) {
            atoms[j] = list[i - j - 1];
        }
        return atoms;
    }

    private final Direction UP = Direction.UP;
    private final Direction DOWN = Direction.DOWN;

}
