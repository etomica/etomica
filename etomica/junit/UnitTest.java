package etomica.junit;

import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesRoot;
import etomica.SpeciesSpheres;
import etomica.SpeciesSpheresMono;
import etomica.SpeciesTree;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.space3d.Space3D;

/**
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Style - Code Templates
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Apr 28, 2005 by kofke
 */
public class UnitTest {

    public static boolean VERBOSE = false;

    /**
     *  
     */
    private UnitTest() {
        super();
        // TODO Auto-generated constructor stub
    }

    public static SpeciesRoot makeStandardSpeciesTree() {
        return makeStandardSpeciesTree(new int[] { 5, 7 }, 3, new int[] { 10,
                10 }, new int[] { 3, 3 }, new int[] { 5, 4, 3 });
    }

    /**
     * Makes tree hierarchy of three species and one or more phases. First
     * species has multiatoms molecules, second species has monatomic molecules
     * (with leaf atoms in molecule layer), and third species has an arbitrary
     * tree structure.
     * 
     * @param n0
     *            number of atoms of species0 in each phase
     * @param nA0
     *            number of atoms per molecule of species0
     * @param n1
     *            number of atoms of species1 in each phase
     * @param n2
     *            number of atoms of species2 in each phase
     * @param n2A
     *            tree specification of species 2
     * @return root of species hierarchy
     */

    public static SpeciesRoot makeStandardSpeciesTree(int[] n0, int nA0,
            int[] n1, int[] n2, int[] n2Tree) {
        Space space = new Space3D();
        Simulation sim = new Simulation(space, new PotentialMaster(space),
                new int[] { 1, 4, 4, 11, 6, 3, 3 });
        Species species0 = null;
        Species species1 = null;
        Species species2 = null;
        int nPhase = 0;
        if (n0 != null) {
            species0 = new SpeciesSpheres(sim, nA0);
            nPhase = n0.length;
        }
        if (n1 != null) {
            species1 = new SpeciesSpheresMono(sim);
            nPhase = n1.length;
        }
        if (n2 != null) {
            species2 = new SpeciesTree(sim, n2Tree);
            nPhase = n2.length;
        }
        for (int i = 0; i < nPhase; i++) {
            Phase phase = new Phase(sim);
            if (species0 != null)
                phase.getAgent(species0).setNMolecules(n0[i]);
            if (species1 != null)
                phase.getAgent(species1).setNMolecules(n1[i]);
            if (species2 != null)
                phase.getAgent(species2).setNMolecules(n2[i]);
        }
        return sim.speciesRoot;
    }

    public static void main(String[] arg) {
        SpeciesRoot root = makeStandardSpeciesTree();
        AtomIteratorTree iterator = new AtomIteratorTree();
        iterator.setRoot(root);
        iterator.setDoAllNodes(true);
        iterator.reset();
        while (iterator.hasNext()) {
            System.out.println(iterator.next().toString());
        }
    }

}