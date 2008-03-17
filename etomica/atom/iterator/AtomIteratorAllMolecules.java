package etomica.atom.iterator;

import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
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
        super(new AtomIteratorArrayListSimple());
    }

    /**
     * Returns a new iterator ready to iterate over the molecules of the given
     * box.
     */
    public AtomIteratorAllMolecules(IBox box) {
        this();
        setBox(box);
    }

    /**
     * Sets the box having the molecules to be returned as iterates.
     */
    public void setBox(IBox box) {
        ((AtomIteratorArrayListSimple)iterator).setList(box.getMoleculeList());
    }
    
    private static final long serialVersionUID = 1L;

    /**
     * main method to test and demonstrate use of this class.
     */
    public static void main(String args[]) {

        Space space = Space2D.getInstance();
        ISimulation sim = new Simulation(space);
        ISpecies species2 = new SpeciesSpheresMono(sim, space);
        ISpecies species1 = new SpeciesSpheres(sim, space, 3);
        ISpecies species0 = new SpeciesSpheres(sim, space, 2);
        sim.getSpeciesManager().addSpecies(species2);
        sim.getSpeciesManager().addSpecies(species1);
        sim.getSpeciesManager().addSpecies(species0);
        IBox box = new Box(sim, space);
        sim.addBox(box);
        box.setNMolecules(species0, 3);
        box.setNMolecules(species1, 2);
        box.setNMolecules(species2, 3);

        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules();

        iterator.setBox(box);
        iterator.reset();
        for (IAtomSet atom = iterator.next(); atom != null; atom = iterator.next()) {
            System.out.println(atom.toString());
        }
        System.out.println();

    }//end of main

}//end of AtomIteratorAllMolecules
