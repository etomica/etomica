package etomica;

/**
 * One-molecule, intra-molecular potential, defining all interatomic potentials within a molecule.
 * All molecules of the associated species (set by speciesIndex) have intramolecular interactions 
 * as defined here.
 */
 
 public abstract class Potential1 implements Simulation.Element, java.io.Serializable {

/**
 * Index of the species that has its intramolecular interactions defined here
 *
 * @see Species#speciesIndex
 */
  int speciesIndex;
  
  private String name;
  private final Simulation parentSimulation;
  private boolean added = false;

  public Potential1(Simulation sim) {
    parentSimulation = sim;
    speciesIndex = 0;
    parentSimulation.register(this);
  }
  
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Potential1.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    
  /**
   * Returns the interatomic potential between the two given atoms.
   * Both atoms should be part of the same molecule
   */
  public abstract Potential getPotential(Atom a1, Atom a2);
  
  /**
   * Accessor method for the speciesIndex, which defines the species that this potential describes
   */
  public final int getSpeciesIndex() {return this.speciesIndex;}  
  /**
   * Accessor method for the speciesIndex, which defines the species that this potential describes
   */
  public final void setSpeciesIndex(int index) {this.speciesIndex = index;}
  
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


