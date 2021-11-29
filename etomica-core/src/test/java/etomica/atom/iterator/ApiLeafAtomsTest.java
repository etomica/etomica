/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.UnitTestUtil;
import etomica.simulation.Simulation;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static etomica.atom.iterator.IteratorTestAbstract.countTest;
import static etomica.atom.iterator.IteratorTestAbstract.testNoIterates;


/**
 * Unit test for ApiLeafAtoms
 *
 * @author David Kofke
 *
 */
public class ApiLeafAtomsTest {

    @Test
    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 0, 6};
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1);
        
        ApiLeafAtoms api = new ApiLeafAtoms();
        
        //test new iterator gives no iterates
        testNoIterates(api);

        int boxCount = sim.getBoxCount();
        for(int i=0; i<boxCount; i++) {
            api.setBox(sim.getBox(i));
            int count = nA0*n0[i] + n1[i];
            count = count*(count-1)/2;
            countTest(api, count);
        }
        
        boolean exceptionThrown = false;
        try {
            api.setBox(null);
        }
        catch (RuntimeException e) {
            exceptionThrown = true;
        }
        Assertions.assertTrue(exceptionThrown);

//        //test documented exceptions
//        boolean exceptionThrown = false;
//        try {
//            new ApiInterspeciesAA(new Species[] {species[0]});
//        } catch(IllegalArgumentException e) {exceptionThrown = true;}
//        assertTrue(exceptionThrown);
//        exceptionThrown = false;
//        try {
//            new ApiInterspeciesAA(new Species[] {species[0], species[0]});
//        } catch(IllegalArgumentException e) {exceptionThrown = true;}
//        assertTrue(exceptionThrown);
//        exceptionThrown = false;
//        try {
//            new ApiInterspeciesAA(new Species[] {species[0], null});
//        } catch(NullPointerException e) {exceptionThrown = true;}
//        assertTrue(exceptionThrown);
//        exceptionThrown = false;
//        try {
//            new ApiInterspeciesAA(null);
//        } catch(NullPointerException e) {exceptionThrown = true;}
//        assertTrue(exceptionThrown);

    }
}
