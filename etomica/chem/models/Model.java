package etomica.chem.models;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.potential.IPotential;
import etomica.simulation.ISimulation;
import etomica.species.Species;

/**
 * Top-level class for a molecular model.
 * @author Andrew Schultz
 */
public abstract class Model implements java.io.Serializable {
	
	/**
     * Returns the species associated with this Model, if it has already been
     * created.  If the species has not been made yet, getSpecies returns null.
	 */
    public Species getSpecies() {
        return species;
    }
    
    /**
     * Creates a species in the given simulation and returns it.
     */
	public final Species makeSpecies(ISimulation sim) {
        if (species == null) {
            species = makeSpeciesInternal(sim);
            sim.getSpeciesManager().addSpecies(species);
            initPotentials(sim);
        }
        return species;
    }

    /**
     * Internal method to be implemented by subclasses to create the actual
     * Species object for the given Simulation.
     */
    protected abstract Species makeSpeciesInternal(ISimulation sim);
	
    /**
     * Internal method to be implemented by subclasses to initialize the
     * intramolecular potentials associated with this model.  Potential
     * objects might be created earlier, but this method gives the subclass
     * an opportunity to create the potentials after the species.  This method
     * will only be called after the species has been created.
     */
    protected abstract void initPotentials(ISimulation sim);
    
    /**
     * Returns an array of objects wrapping bonding Potentials and the 
     * AtomsetIteratorBasisDependents that return atoms appropriate for those
     * potentials.
     */
	public abstract PotentialAndIterator[] getPotentials();
    
    private Species species;
    protected PotentialAndIterator[] potentialsAndIterators;

    /**
     * Wrapper class for a Potential and an AtomsetIteratorBasisDependent that
     * is associated with the potential.
     */
    public static class PotentialAndIterator {
        protected PotentialAndIterator(IPotential potential, AtomsetIteratorBasisDependent iterator) {
            this.potential = potential;
            this.iterator = iterator;
        }
        
        public IPotential getPotential() {
            return potential;
        }
        
        public AtomsetIteratorBasisDependent getIterator() {
            return iterator;
        }
        
        private final IPotential potential;
        private final AtomsetIteratorBasisDependent iterator;
    }
}
