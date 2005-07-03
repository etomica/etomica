package etomica.junit;

import etomica.AtomFactory;
import etomica.AtomTypeGroup;
import etomica.AtomTypeLeaf;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesRoot;
import etomica.SpeciesSpheres;
import etomica.SpeciesSpheresMono;
import etomica.SpeciesTree;
import etomica.atom.AtomFactoryHetero;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.space3d.Space3D;

/**
 * Contains some convenience methods and fields useful for implementing unit
 * tests.
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
     * Private to prevent instantiation
     */
    private UnitTest() {
        super();
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
     *            tree specification of species 2, e.g., {2,4} indicates that
     *            each molecule has 2 subgroups, each with 4 atoms (such as
     *            CH3-CH3)
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

    /**
     * Makes tree hierarchy of one or more species in a single phase. Number of
     * species if determined by length of nMolecules array. Each molecule is
     * heterogeneous, formed from atoms of different types. For example:
     * <ul>
     * <li>nMolecules = {5,3} will form two species, with five molecules of
     * speciesA and 3 molecule of SpeciesB.
     * <li>Then with nAtoms = {{2,1,4},{2,3}}, a speciesA molecule will contain
     * 2+1+4 = 7 atoms, with 2 of type-a, 1 of type-b, and 4 of type-c, and a
     * speciesB molecule will contain 2+3 = 5 atoms, with 2 of type-d and 3 of
     * type-e. All types are different instance of AtomTypeSphere.
     * </ul>
     * 
     * @param nMolecules
     *            number of molecules made of each species. Number of species is
     *            determined by length of this array.
     * @param nAtoms
     *            first index corresponds to species, so nAtoms.length should
     *            equal nMolecules.length. Number of types in a species-j
     *            molecule is given by the length of the nAtoms[j] subarray, and
     *            the elements of this array give the number of atoms of each
     *            type used to form a molecule.
     * @return root of the species hierarchy
     */
    public static SpeciesRoot makeMultitypeSpeciesTree(int[] nMolecules,
            int[][] nAtoms) {
        Space space = new Space3D();
        Simulation sim = new Simulation(space, new PotentialMaster(space),
                new int[] { 1, 4, 4, 11, 6, 3, 3 });
        //        new SpeciesSpheres(sim);
        for (int i = 0; i < nMolecules.length; i++) {
            AtomTypeGroup agentType = Species.makeAgentType(sim);
            AtomFactoryHetero factory = new AtomFactoryHetero(space,
                    AtomSequencerFactory.SIMPLE, agentType);
            AtomFactory[] childFactories = new AtomFactory[nAtoms[i].length];
            for (int j = 0; j < childFactories.length; j++) {
                AtomTypeLeaf atomType = new AtomTypeSphere(
                        (AtomTypeGroup) factory.getType());
                childFactories[j] = new AtomFactoryMono(space, atomType,
                        AtomSequencerFactory.SIMPLE);
            }
            factory.setChildFactory(childFactories, nAtoms[i]);
            Species species = new Species(sim, factory, agentType);
            species.setNMolecules(nMolecules[i]);
        }
        new Phase(sim);
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