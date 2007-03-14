package etomica.junit.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.ApiIntraspecies1A;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.species.Species;


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
public class ApiIntraspecies1ATest extends IteratorTestAbstract {

    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 1, 6};
        int[] n2 = new int[] {1, 7, 2};
        int[] n2Tree = new int[] {3,4};
        SpeciesRoot root = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2, n2Tree);
        
        Species[] species = new Species[3];
        species[0] = root.getDescendant(new int[] {0,0}).getType().getSpecies();
        species[1] = root.getDescendant(new int[] {0,1}).getType().getSpecies();
        species[2] = root.getDescendant(new int[] {0,2}).getType().getSpecies();

        phaseTest(root, species, 0);
        phaseTest(root, species, 1);
        phaseTest(root, species, 2);
        
        ApiIntraspecies1A api = new ApiIntraspecies1A(new Species[] {species[0], species[0]});
        
        //test new iterator gives no iterates
        testNoIterates(api);

        //one species has no molecules
        api.setPhase(root.getDescendant(new int[] {2}).parentPhase());
        api.setTarget(root.getDescendant(new int[] {2,1,3}));
        testNoIterates(api);
        //target not one of species
        api = new ApiIntraspecies1A(new Species[] {species[1],species[1]});
        api.setPhase(root.getDescendant(new int[] {0}).parentPhase());
        api.setTarget(root.getDescendant(new int[] {0,0,3}));
        testNoIterates(api);
        //target one of species but in different phase
        api = new ApiIntraspecies1A(new Species[] {species[1],species[1]});
        api.setPhase(root.getDescendant(new int[] {0}).parentPhase());
        api.setTarget(root.getDescendant(new int[] {1,1,0}));
        testNoIterates(api);
        
        //test documented exceptions
        Atom target = null;
        boolean exceptionThrown = false;
        try {
            api.setTarget(target);
        } catch(NullPointerException e) {exceptionThrown = true;}
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
    private void phaseTest(SpeciesRoot root, Species[] species, int phaseIndex) {
        speciesTestForward(root, species, phaseIndex, 0);
        speciesTestForward(root, species, phaseIndex, 1);
        speciesTestForward(root, species, phaseIndex, 2);
    }

    /**
     * Test iteration in various directions with different targets.  Iterator constructed with
     * index of first species less than index of second.
     */
    private void speciesTestForward(SpeciesRoot root, Species[] species, int phaseIndex, int species0Index) {

        ApiIntraspecies1A api = new ApiIntraspecies1A(new Species[] {species[species0Index], species[species0Index]});
        Phase phase = root.getDescendant(new int[] {phaseIndex}).parentPhase();
        AtomsetAction speciesTest = new SpeciesTestAction(species[species0Index], species[species0Index]);
        Atom target = null;
        Atom targetMolecule = null;
        //test no iterates if no target
        api.setPhase(phase);
        Atom[] molecules0 = phase.getAgent(species[species0Index]).getChildList().toArray();
        int[] nMolecules = new int[] {molecules0.length};
        testNoIterates(api);
        
        if(nMolecules[0] == 0) return;

        //species0 target; any direction
        
        target = root.getDescendant(new int[] {phaseIndex,species0Index,nMolecules[0]/2});
        targetMolecule = target;
        api.setTarget(target);
        testApiIterates(api,targetMolecule,upMolecules(targetMolecule,molecules0), dnMolecules(targetMolecule, molecules0));
        api.allAtoms(speciesTest);

        //species0 target; up
        target = root.getDescendant(new int[] {phaseIndex,species0Index,nMolecules[0]/2});
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(UP);
        testApiIterates(api,UP, targetMolecule,upMolecules(targetMolecule,molecules0));
        api.allAtoms(speciesTest);

        //species0 target; down
        target = root.getDescendant(new int[] {phaseIndex,species0Index,nMolecules[0]/2});
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(DOWN);
        testApiIteratesSwap(api,targetMolecule,dnMolecules(targetMolecule,molecules0));
        
        //species0 leafAtom target; any direction
        if(species0Index != 1) {
            target = root.getDescendant(new int[] {phaseIndex,species0Index,nMolecules[0]/2,1});
            targetMolecule = target.parentGroup();
            api.setTarget(target);
            api.setDirection(UP);
            testApiIterates(api,UP, targetMolecule,upMolecules(targetMolecule, molecules0));
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
        public void actionPerformed(AtomSet atomSet) {
            assertTrue(atomSet.getAtom(0).getType().getSpecies() == species0);
            assertTrue(atomSet.getAtom(1).getType().getSpecies() == species1);
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
    
    private final IteratorDirective.Direction UP = IteratorDirective.Direction.UP;
    private final IteratorDirective.Direction DOWN = IteratorDirective.Direction.DOWN;

}
