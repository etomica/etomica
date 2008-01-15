package etomica.junit.atom.iterator;

import etomica.atom.iterator.ApiLeafAtoms;
import etomica.junit.UnitTestUtil;
import etomica.box.Box;
import etomica.simulation.ISimulation;


/**
 * Unit test for ApiLeafAtoms
 *
 * @author David Kofke
 *
 */
public class ApiLeafAtomsTest extends IteratorTestAbstract {

    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 0, 6};
        ISimulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1);
        
        ApiLeafAtoms api = new ApiLeafAtoms();
        
        //test new iterator gives no iterates
        testNoIterates(api);
        
        Box[] box = sim.getBoxs();
        
        for(int i=0; i<box.length; i++) {
            api.setBox(box[i]);
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
        assertTrue(exceptionThrown);

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
