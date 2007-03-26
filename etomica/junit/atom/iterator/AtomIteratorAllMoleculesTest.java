package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
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
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2,
                n2Tree);

        Species[] species = sim.getSpeciesManager().getSpecies();

        for(int i=0; i<n0.length; i++) {
            phaseTest(sim.getPhases()[i], species);
        }

    }

    /**
     * Performs tests on different species combinations in a particular phase.
     */
    private void phaseTest(Phase phase, Species[] species) {
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules();

        iterator.setPhase(phase);
        
        AtomArrayList moleculeList = new AtomArrayList();
        for(int i=0; i<species.length; i++) {
            AtomArrayList molecules = phase.getAgent(species[i]).getChildList();
            for (int j=0; j<molecules.size(); j++) {
                moleculeList.add(molecules.get(j));
            }
        }
        
        LinkedList list = testIterates(iterator, moleculeList.toArray());
        assertEquals(list.size(), phase.moleculeCount());
    }
}
