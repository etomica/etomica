package etomica.junit.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetActionAdapter;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.iterator.ApiIntraspecies1A;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;


/**
 * Unit test for Apiintraspecies1A
 *
 * @author David Kofke
 *
 */
public class ApiIntraspecies1ATest extends IteratorTestAbstract {

    public void testIterator() {
        
        int[] n0 = new int[] {10, 1, 0};
        int nA0 = 5;
        int[] n1 = new int[] {5, 1, 6};
        int[] n2 = new int[] {1, 7, 2};
        int[] n2Tree = new int[] {3,4};
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(n0, nA0, n1, n2, n2Tree);
        
        Species[] species = sim.getSpeciesManager().getSpecies();

        phaseTest(sim.getPhases()[0], species);
        phaseTest(sim.getPhases()[1], species);
        phaseTest(sim.getPhases()[2], species);
        
        ApiIntraspecies1A api = new ApiIntraspecies1A(new Species[] {species[0], species[0]});
        
        //test new iterator gives no iterates
        testNoIterates(api);

        //test documented exceptions
        IAtom target = null;
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
    private void phaseTest(Phase phase, Species[] species) {
        speciesTestForward(phase, species[0]);
        speciesTestForward(phase, species[1]);
        speciesTestForward(phase, species[2]);
    }

    /**
     * Test iteration in various directions with different targets.  Iterator constructed with
     * index of first species less than index of second.
     */
    private void speciesTestForward(Phase phase, Species species) {

        ApiIntraspecies1A api = new ApiIntraspecies1A(new Species[] {species, species});
        AtomsetAction speciesTest = new SpeciesTestAction(species, species);
        IAtom target = null;
        IAtom targetMolecule = null;
        //test no iterates if no target
        api.setPhase(phase);
        IAtom[] molecules0 = ((AtomArrayList)phase.getAgent(species).getChildList()).toArray();
        int[] nMolecules = new int[] {molecules0.length};
        testNoIterates(api);
        
        if(nMolecules[0] == 0) return;

        //species0 target; any direction
        
        target = phase.getAgent(species).getChildList().getAtom(nMolecules[0]/2);
        targetMolecule = target;
        api.setTarget(target);
        testApiIterates(api,targetMolecule,upMolecules(targetMolecule,molecules0), dnMolecules(targetMolecule, molecules0));
        api.allAtoms(speciesTest);

        //species0 target; up
        target = phase.getAgent(species).getChildList().getAtom(nMolecules[0]/2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(UP);
        testApiIterates(api,UP, targetMolecule,upMolecules(targetMolecule,molecules0));
        api.allAtoms(speciesTest);

        //species0 target; down
        target = phase.getAgent(species).getChildList().getAtom(nMolecules[0]/2);
        targetMolecule = target;
        api.setTarget(target);
        api.setDirection(DOWN);
        testApiIteratesSwap(api,targetMolecule,dnMolecules(targetMolecule,molecules0));
        
        //species0 leafAtom target; any direction
        if(species instanceof SpeciesSpheres) {
            target = ((IAtomGroup)phase.getAgent(species).getChildList().getAtom(nMolecules[0]/2)).getChildList().getAtom(1);
            targetMolecule = target.getParentGroup();
            api.setTarget(target);
            api.setDirection(UP);
            testApiIterates(api,UP, targetMolecule,upMolecules(targetMolecule, molecules0));
            api.allAtoms(speciesTest);
        }

        boolean exceptionThrown = false;
        try {
            api.setPhase(null);
        }
        catch (RuntimeException e) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
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
    
    private IAtom[] upMolecules(IAtom target, IAtom[] list) {
        int i;
        for(i=0; i<list.length; i++) {
            if(list[i] == target) break;
        }
        IAtom[] atoms = new IAtom[list.length-i-1];
        for(int j=0; j<atoms.length; j++) {
            atoms[j] = list[i+j+1];
        }
        return atoms;
    }
    private IAtom[] dnMolecules(IAtom target, IAtom[] list) {
        int i;
        for(i=0; i<list.length; i++) {
            if(list[i] == target) break;
        }
        IAtom[] atoms = new IAtom[i];
        for(int j=0; j<i; j++) {
            atoms[j] = list[i-j-1];
        }
        return atoms;
    }
    
    private final IteratorDirective.Direction UP = IteratorDirective.Direction.UP;
    private final IteratorDirective.Direction DOWN = IteratorDirective.Direction.DOWN;

}
