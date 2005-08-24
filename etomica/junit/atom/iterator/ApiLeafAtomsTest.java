package etomica.junit.atom.iterator;

import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.junit.UnitTest;
import etomica.phase.Phase;


/**
 * Unit test for ApiLeafAtoms
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 28, 2005 by kofke
 */
public class ApiLeafAtomsTest extends IteratorTest {

    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 0, 6};
        int[] n2 = new int[] {1, 7, 2};
        int[] n2Tree = new int[] {3,4};
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(n0, nA0, n1, n2, n2Tree);
        AtomTreeNodeGroup rootNode = (AtomTreeNodeGroup)root.node;
        
        ApiLeafAtoms api = new ApiLeafAtoms();
        
        //test new iterator gives no iterates
        testNoIterates(api);
        
        Phase[] phase = new Phase[3];
        phase[0] = rootNode.getDescendant(new int[] {0}).node.parentPhase();
        phase[1] = rootNode.getDescendant(new int[] {1}).node.parentPhase();
        phase[2] = rootNode.getDescendant(new int[] {2}).node.parentPhase();
        
        for(int i=0; i<phase.length; i++) {
            api.setPhase(phase[i]);
            int count = nA0*n0[i] + n1[i] + n2[i]*n2Tree[0]*n2Tree[1];
            count = count*(count-1)/2;
            countTest(api, count);
        }
        
        api.setPhase(null);
        testNoIterates(api);


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
