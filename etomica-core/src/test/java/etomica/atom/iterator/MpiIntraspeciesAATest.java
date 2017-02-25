/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.UnitTestUtil;
import etomica.atom.MoleculesetAction;


/**
 * Unit test for ApiIntraspeciesAA
 *
 * @author David Kofke
 *
 */

public class MpiIntraspeciesAATest extends MoleculeIteratorTestAbstract {

    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 1, 6};
        ISimulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1);
        
        ISpecies[] species = new ISpecies[sim.getSpeciesCount()];
        for(int i = 0; i < sim.getSpeciesCount(); i++) {
        	species[i] = sim.getSpecies(i);
        }
        boxTest(sim.getBox(0), species);
        boxTest(sim.getBox(1), species);
        boxTest(sim.getBox(2), species);
        
        
        //test new iterator gives no iterates
        MpiIntraspeciesAA api = new MpiIntraspeciesAA(species[0]);

        //test documented exceptions
        boolean exceptionThrown = false;
        try {
            // using iterator without calling setBox
            testNoIterates(api);
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        try {
            // null box
            api.setBox(null);
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        try {
            // null Species
            new MpiIntraspeciesAA((ISpecies)null);
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
    }
    
    /**
     * Performs tests on different species combinations in a particular box.
     */
    private void boxTest(IBox box, ISpecies[] species) {
        speciesTestForward(box, species, 0);
        speciesTestForward(box, species, 1);
    }

    /**
     * Test iteration in various directions with different targets.  Iterator constructed with
     * index of first species less than index of second.
     */
    private void speciesTestForward(IBox box, ISpecies[] species, int species0Index) {
        MpiIntraspeciesAA api = new MpiIntraspeciesAA(species[species0Index]);
        MoleculesetAction speciesTest = new SpeciesTestAction(species[species0Index], species[species0Index]);

        api.setBox(box);
        int nMolecules = box.getMoleculeList(species[species0Index]).getMoleculeCount();
        
        int count = nMolecules * (nMolecules - 1) / 2;
        
        countTest(api, count);
        MoleculeIteratorTestAbstract.allAtoms(api,speciesTest);

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

    private class SpeciesTestAction implements MoleculesetAction {
        final ISpecies species0, species1;
        public SpeciesTestAction(ISpecies species0, ISpecies species1) {
            this.species0 = species0;
            this.species1 = species1;
        }
        public void actionPerformed(IMoleculeList atomSet) {
            assertTrue(atomSet.getMolecule(0).getType() == species0);
            assertTrue(atomSet.getMolecule(1).getType() == species1);
        }
    }
    
}
