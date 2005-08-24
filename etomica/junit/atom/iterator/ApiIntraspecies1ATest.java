package etomica.junit.atom.iterator;

import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Species;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomsetArray;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.ApiIntraspecies1A;
import etomica.junit.UnitTest;


/**
 * Unit test for Apiintraspecies1A
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 28, 2005 by kofke
 */
public class ApiIntraspecies1ATest extends IteratorTest {

    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 1, 6};
        int[] n2 = new int[] {1, 7, 2};
        int[] n2Tree = new int[] {3,4};
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(n0, nA0, n1, n2, n2Tree);
        AtomTreeNodeGroup rootNode = (AtomTreeNodeGroup)root.node;
        
        Species[] species = new Species[3];
        species[0] = rootNode.getDescendant(new int[] {0,0}).type.getSpecies();
        species[1] = rootNode.getDescendant(new int[] {0,1}).type.getSpecies();
        species[2] = rootNode.getDescendant(new int[] {0,2}).type.getSpecies();

        phaseTest(rootNode, species, 0);
        phaseTest(rootNode, species, 1);
        phaseTest(rootNode, species, 2);
        
        ApiIntraspecies1A api = new ApiIntraspecies1A(new Species[] {species[0], species[0]});
        
        //test new iterator gives no iterates
        testNoIterates(api);

        //one species has no molecules
        api.setPhase(rootNode.getDescendant(new int[] {2}).node.parentPhase());
        api.setTarget(rootNode.getDescendant(new int[] {2,1,3}));
        testNoIterates(api);
        //target not one of species
        api = new ApiIntraspecies1A(new Species[] {species[1],species[1]});
        api.setPhase(rootNode.getDescendant(new int[] {0}).node.parentPhase());
        api.setTarget(rootNode.getDescendant(new int[] {0,0,3}));
        testNoIterates(api);
        //target one of species but in different phase
        api = new ApiIntraspecies1A(new Species[] {species[1],species[1]});
        api.setPhase(rootNode.getDescendant(new int[] {0}).node.parentPhase());
        api.setTarget(rootNode.getDescendant(new int[] {1,1,0}));
        testNoIterates(api);
        //target atom is null
        api = new ApiIntraspecies1A(new Species[] {species[0],species[0]});
        api.setPhase(rootNode.getDescendant(new int[] {0}).node.parentPhase());
        api.setTarget(new AtomsetArray(1));
        testNoIterates(api);
        
        
        //test documented exceptions
        AtomSet target = AtomSet.NULL;
        boolean exceptionThrown = false;
        try {
            api.setTarget(target);
        } catch(IllegalArgumentException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            api.setTarget(null);
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiIntraspecies1A(new Species[] {species[0]});
        } catch(IllegalArgumentException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiIntraspecies1A(new Species[] {species[0], species[1]});
        } catch(IllegalArgumentException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiIntraspecies1A(new Species[] {species[0], null});
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        try {
            new ApiIntraspecies1A((Species)null);
        } catch(NullPointerException e) {exceptionThrown = true;}
        try {
            new ApiIntraspecies1A((Species[])null);
        } catch(NullPointerException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);

        
    }
    
    /**
     * Performs tests on different species combinations in a particular phase.
     */
    private void phaseTest(AtomTreeNodeGroup rootNode, Species[] species, int phaseIndex) {
        speciesTestForward(rootNode, species, phaseIndex, 0);
        speciesTestForward(rootNode, species, phaseIndex, 1);
        speciesTestForward(rootNode, species, phaseIndex, 2);
    }

    /**
     * Test iteration in various directions with different targets.  Iterator constructed with
     * index of first species less than index of second.
     */
    private void speciesTestForward(AtomTreeNodeGroup rootNode, Species[] species, int phaseIndex, int species0Index) {

        ApiIntraspecies1A api = new ApiIntraspecies1A(new Species[] {species[species0Index], species[species0Index]});
        Phase phase = rootNode.getDescendant(new int[] {phaseIndex}).node.parentPhase();
        AtomsetAction speciesTest = new SpeciesTestAction(species[species0Index], species[species0Index]);
        AtomSet target = AtomSet.NULL;
        Atom targetMolecule = null;
        //test no iterates if no target
        api.setPhase(phase);
        Atom[] molecules0 = ((AtomTreeNodeGroup)phase.getAgent(species[species0Index]).node).childList.toArray();
        int[] nMolecules = new int[] {molecules0.length};
        testNoIterates(api);
        
        if(nMolecules[0] == 0) return;

        //species0 target; any direction
        
        target = rootNode.getDescendant(new int[] {phaseIndex,species0Index,nMolecules[0]/2});
        targetMolecule = (Atom)target;
        api.setTarget(target);
        testApiIterates(api,targetMolecule,upMolecules(targetMolecule,molecules0), dnMolecules(targetMolecule, molecules0));
        api.allAtoms(speciesTest);

        //species0 target; up
        target = rootNode.getDescendant(new int[] {phaseIndex,species0Index,nMolecules[0]/2});
        targetMolecule = (Atom)target;
        api.setTarget(target);
        api.setDirection(UP);
        testApiIterates(api,targetMolecule,upMolecules(targetMolecule,molecules0));
        api.allAtoms(speciesTest);

        //species0 target; down
        target = rootNode.getDescendant(new int[] {phaseIndex,species0Index,nMolecules[0]/2});
        targetMolecule = (Atom)target;
        api.setTarget(target);
        api.setDirection(DOWN);
        testApiIteratesSwap(api,targetMolecule,dnMolecules(targetMolecule,molecules0));
        
        //species0 leafAtom target; any direction
        if(species0Index != 1) {
            target = rootNode.getDescendant(new int[] {phaseIndex,species0Index,nMolecules[0]/2,1});
            targetMolecule = ((Atom)target).node.parentGroup();
            api.setTarget(target);
            api.setDirection(UP);
            testApiIterates(api,targetMolecule,upMolecules(targetMolecule, molecules0));
            api.allAtoms(speciesTest);
        }

        api.setPhase(null);
        testNoIterates(api);
    }

    private class SpeciesTestAction extends AtomsetActionAdapter {
        final Species species0, species1;
        public SpeciesTestAction(Species species0, Species species1) {
            this.species0 = species0;
            this.species1 = species1;
        }
        public void actionPerformed(AtomSet atoms) {
            assertTrue(atoms.getAtom(0).type.getSpecies() == species0);
            assertTrue(atoms.getAtom(1).type.getSpecies() == species1);
        }
    }
    
    private Atom[] upMolecules(Atom target, Atom[] list) {
        int i;
        for(i=0; i<list.length; i++) {
            if(list[i] == target) break;
        }
        Atom[] atoms = new Atom[list.length-i-1];
        for(int j=0; j<atoms.length; j++) {
            atoms[j] = list[i+j+1];
        }
        return atoms;
    }
    private Atom[] dnMolecules(Atom target, Atom[] list) {
        int i;
        for(i=0; i<list.length; i++) {
            if(list[i] == target) break;
        }
        Atom[] atoms = new Atom[i];
        for(int j=0; j<i; j++) {
            atoms[j] = list[i-j-1];
        }
        return atoms;
    }
    
    private final IteratorDirective.Direction UP = IteratorDirective.UP;
    private final IteratorDirective.Direction DOWN = IteratorDirective.DOWN;

}
