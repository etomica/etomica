package etomica.atom.iterator;

import etomica.atom.AtomSet;
import etomica.box.Box;
import etomica.simulation.ISimulation;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;

/**
 * Iterator for all the molecules in a box. Loops over all those atoms that
 * lie just below the species agents in the atom tree hierarchy. To iterate the
 * molecules of just one species, use AtomIteratorMolecule.
 * 
 * @author David Kofke
 * @since 02.02.16
 */

public class AtomIteratorAllMolecules extends AtomIteratorAdapter 
            implements AtomIteratorBoxDependent {

    public AtomIteratorAllMolecules() {
        super(new AtomIteratorTreeBox(2));
    }

    /**
     * Returns a new iterator ready to iterate over the molecules of the given
     * box.
     */
    public AtomIteratorAllMolecules(Box box) {
        this();
        setBox(box);
    }

    /**
     * Sets the box having the molecules to be returned as iterates.
     */
    public void setBox(Box box) {
        ((AtomIteratorBoxDependent)iterator).setBox(box);
    }
    
    private static final long serialVersionUID = 1L;

    /**
     * main method to test and demonstrate use of this class.
     */
    public static void main(String args[]) {

        ISimulation sim = new Simulation();
        Species species2 = new SpeciesSpheresMono(sim);
        Species species1 = new SpeciesSpheres(sim, 3);
        Species species0 = new SpeciesSpheres(sim, 2);
        sim.getSpeciesManager().addSpecies(species2);
        sim.getSpeciesManager().addSpecies(species1);
        sim.getSpeciesManager().addSpecies(species0);
        Box box = new Box(sim);
        sim.addBox(box);
        box.setNMolecules(species0, 3);
        box.setNMolecules(species1, 2);
        box.setNMolecules(species2, 3);

        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules();

        iterator.setBox(box);
        iterator.reset();
        for (AtomSet atom = iterator.next(); atom != null; atom = iterator.next()) {
            System.out.println(atom.toString());
        }
        System.out.println();

    }//end of main

}//end of AtomIteratorAllMolecules
