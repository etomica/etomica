package etomica;

/**
 * Intermolecular pair potential formed from a single interatomic potential
 * All atoms of two molecules interact according to the single potential used
 * to define this intermolecular potential
 */
public class P2SimpleWrapper extends Potential2 {

  private Potential onlyPotential;

  public P2SimpleWrapper(Potential p) {
    this(Simulation.instance, p);
  }
  public P2SimpleWrapper(Simulation sim, Potential p) {
    super(sim);
    onlyPotential = p;
  }
  
  /**
   * Implementation of the Potential2 interface to return the atomic potential 
   * that governs the interaction between the given atoms.
   */
  public final Potential getPotential(Atom a1, Atom a2) {return onlyPotential;}
  
  /**
   * Accessor method for the interatomic potential
   */
  public final Potential getOnlyPotential() {return onlyPotential;}
  /**
   * Accessor method for the interatomic potential
   */
  public final void setOnlyPotential(Potential p) {onlyPotential = p;}

}


