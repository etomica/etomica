package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.junit.UnitTestUtil;
import etomica.box.Box;
import etomica.simulation.ISimulation;
import etomica.species.Species;

/**
 * Unit test for ApiIntraspeciesAA
 * 
 * @author David Kofke
 *  
 */
public class AtomIteratorAllMoleculesTest extends IteratorTestAbstract {

    public void testIterator() {

        int[] n0 = new int[] { 10, 1, 0, 0, 0};
        int nA0 = 5;
        int[] n1 = new int[] { 5, 1, 6, 0, 1};
        int[] n2 = new int[] { 1, 7, 2, 0, 0};
        int[] n2Tree = new int[] { 3, 4 };
        ISimulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2,
                n2Tree);

        Species[] species = sim.getSpeciesManager().getSpecies();

        for(int i=0; i<n0.length; i++) {
            boxTest(sim.getBoxs()[i], species);
        }

    }

    /**
     * Performs tests on different species combinations in a particular box.
     */
    private void boxTest(Box box, Species[] species) {
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules();

        iterator.setBox(box);
        
        AtomArrayList moleculeList = new AtomArrayList();
        for(int i=0; i<species.length; i++) {
            AtomSet molecules = box.getMoleculeList(species[i]);
            for (int j=0; j<molecules.getAtomCount(); j++) {
                moleculeList.add(molecules.getAtom(j));
            }
        }
        
        LinkedList list = testIterates(iterator, moleculeList.toArray());
        assertEquals(list.size(), box.moleculeCount());
    }
}
