/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import java.util.LinkedList;

import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.simulation.Simulation;
import etomica.api.ISpecies;
import etomica.atom.MoleculeArrayList;
import etomica.UnitTestUtil;

/**
 * Unit test for ApiIntraspeciesAA
 * 
 * @author David Kofke
 *  
 */
public class MoleculeIteratorAllMoleculesTest extends MoleculeIteratorTestAbstract {

    public void testIterator() {

        int[] n0 = new int[] { 10, 1, 0, 0, 0};
        int nA0 = 5;
        int[] n1 = new int[] { 5, 1, 6, 0, 1};
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1);

        ISpecies[] species = new ISpecies[sim.getSpeciesCount()];
        for(int i = 0; i < sim.getSpeciesCount(); i++) {
        	species[i] = sim.getSpecies(i);
        }
        for(int i=0; i<n0.length; i++) {
            boxTest(sim.getBox(i), species);
        }

    }

    /**
     * Performs tests on different species combinations in a particular box.
     */
    private void boxTest(Box box, ISpecies[] species) {
        MoleculeIteratorAllMolecules iterator = new MoleculeIteratorAllMolecules();

        iterator.setBox(box);
        
        MoleculeArrayList moleculeList = new MoleculeArrayList();
        for(int i=0; i<species.length; i++) {
            IMoleculeList molecules = box.getMoleculeList(species[i]);
            for (int j=0; j<molecules.getMoleculeCount(); j++) {
                moleculeList.add(molecules.getMolecule(j));
            }
        }
        
        LinkedList list = testIterates(iterator, moleculeList.toMoleculeArray());
        assertEquals(list.size(), box.getMoleculeList().getMoleculeCount());
    }
}
