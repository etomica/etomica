package etomica.atom.iterator;

import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;

/**
 * Iterator for all the molecules in a phase. Loops over all those atoms that
 * lie just below the species agents in the atom tree hierarchy. To iterate the
 * molecules of just one species, use AtomIteratorMolecule.
 * 
 * @author David Kofke
 * @since 02.02.16
 */

public class AtomIteratorAllMolecules extends AtomIteratorAdapter 
            implements AtomIteratorPhaseDependent {

    public AtomIteratorAllMolecules() {
        super(new AtomIteratorTree(2));
        treeIterator = (AtomIteratorTree) iterator;
        setPhase(null);
    }

    /**
     * Returns a new iterator ready to iterate over the molecules of the given
     * phase.
     */
    public AtomIteratorAllMolecules(Phase phase) {
        this();
        setPhase(phase);
    }

    /**
     * Sets the phase having the molecules to be returned as iterates.
     */
    public void setPhase(Phase phase) {
        treeIterator.setRoot((phase == null) ? null : phase.getSpeciesMaster());
    }

    private static final long serialVersionUID = 1L;
    private final AtomIteratorTree treeIterator;

    /**
     * main method to test and demonstrate use of this class.
     */
    public static void main(String args[]) {

        Simulation sim = new Simulation();
        Species species2 = new SpeciesSpheresMono(sim);
        Species species1 = new SpeciesSpheres(sim, 3);
        Species species0 = new SpeciesSpheres(sim, 2);
        sim.getSpeciesRoot().addSpecies(species2);
        sim.getSpeciesRoot().addSpecies(species1);
        sim.getSpeciesRoot().addSpecies(species0);
        Phase phase = new Phase(sim);
        phase.getAgent(species0).setNMolecules(3);
        phase.getAgent(species1).setNMolecules(2);
        phase.getAgent(species2).setNMolecules(3);

        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules();
        System.out.println(iterator.hasNext());

        iterator.setPhase(phase);
        System.out.println(iterator.hasNext());
        iterator.reset();
        while (iterator.hasNext())
            System.out.println(iterator.next().toString());
        System.out.println();

    }//end of main

}//end of AtomIteratorAllMolecules
