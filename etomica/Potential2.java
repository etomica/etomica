package etomica; 

/**
 * Two-molecule, inter-molecular potential, defining all interatomic potentials of two different molecules.
 * All molecules of the two associated species (set by speciesIndex) have intermolecular interactions 
 * as defined here.
 */
public abstract class Potential2 implements Simulation.Element, java.io.Serializable {

/**
 * Index of one of the species
 *
 * @see Species#speciesIndex
 */
  int species1Index;
/**
 * Index of one of the species
 *
 * @see Species#speciesIndex
 */ 
  int species2Index;

  private String name;
  private final Simulation parentSimulation;
  private boolean added = false;

/**
 * Range of potential; molecules separated by more than this amount do not interact
 * Enforcement of cutoff is responsibility of subclasses
 */
  double potentialCutoff;

  public Potential2(Simulation sim) {
    parentSimulation = sim;
    species1Index = species2Index = 0;
    parentSimulation.register(this);
  }

    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Potential2.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    
 /**
  * Returns the inter-atomic potential for the two given atoms
  * In many cases, the same potential governs interactions of all atoms of two molecules.
  */
  public abstract Potential getPotential(Atom a1, Atom a2);
 
 /**
  * Accessor method for potentialCutoff
  * @see #potentialCutoff
  */
  public final void setPotentialCutoff(double d) {
    potentialCutoff = d;
  }
  
 /**
  * Accessor method for species index
  * 
  * @see Species#speciesIndex
  */
  public int getSpecies1Index() {return this.species1Index;}  
 /**
  * Accessor method for species index
  * 
  * @see Species#speciesIndex
  */
  public void setSpecies1Index(int index) {this.species1Index = index;}
  
 /**
  * Accessor method for species index
  * 
  * @see Species#speciesIndex
  */
  public int getSpecies2Index() {return this.species2Index;}
 /**
  * Accessor method for species index
  * 
  * @see Species#speciesIndex
  */
  public void setSpecies2Index(int index) {this.species2Index = index;}
  
  /**
   * Accessor method for both species indices
   */
  public void setSpeciesIndex(int i1, int i2) {
    species1Index = i1;
    species2Index = i2;
  }
  
    /**
     * Accessor method of the name of this object
     * 
     * @return The given name
     */
    public final String getName() {return name;}

    /**
     * Method to set the name of this object
     * 
     * @param name The name string to be associated with this object
     */
    public final void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the object
     */
    public String toString() {return getName();}  //override Object method
          
}


