package simulate;
import java.io.*;
import java.util.*;
import java.awt.Graphics;

/**
 * A Molecule is a collection of Atoms.  Molecules differ in the number
 * and (rarely) type of Atoms that form them.  The number of Atoms in a Molecule
 * is set via an instance variable; Molecules composed of specialized 
 * Atom types (formed as a subclass of Atom) must be defined as a subclass 
 * of Molecule.
 * All Molecules in a Phase are ordered in a linked list.
 * Many instances of Molecules of the same type are collected together
 * in a Species.
 *
 * @author David Kofke
 * @author C. Daniel Barnes
 * @see Atom
 * @see Species
 */

public class Molecule implements Serializable {

  /**
   * The total mass of the molecule; the sum of all atom masses (amu)
   */
  double mass;
  
  /**
   * Collection of molecules used as a neighbor list.  Use is strongly
   * discouraged as present neighbor-list structure is likely to change
   * in future version
   */
  Vector neighbors;
  
  private final double[] r0 = new double[Space.D];
  private final double[] r = new double[Space.D];
  private final double[] dr = new double[Space.D];
  final double[] p = new double[Space.D];  //maybe delete p and f some day
  final double[] f = new double[Space.D];
  
 /**
  * Instance of the species in which this molecule resides
  * Assigned in the Molecule constructor.
  * @see Species#makeMolecules
  */
  Species parentSpecies;
  
 /**
  * Next molecule in linked list of molecules
  * @see #setNextMolecule
  * @see Species#orderMolecules
  */
  Molecule nextMolecule;
 
 /**
  * Previous molecule in linked list of molecules
  * @set #setNextMolecule
  */
  Molecule previousMolecule;
  
 /**
  * Atoms in this molecule are part of a linked list of all atoms
  * in the simulation.  This is the first such atom in this molecule.
  */
  Atom firstAtom;
  
 /**
  * Last atom in the molecule.
  * @see firstAtom
  */
  Atom lastAtom;
  
 /**
  * Array of atoms in this molecule.  This data structure may be
  * eliminated in future revisions in lieu of just the linked list.
  */
//  Atom[] atom;
  
 /**
  * Number of atoms in this molecule.  Assigned by species when invoking
  * the molecule's constructor
  * @see setNAtoms
  * @see Species#makeMolecules
  */
  int nAtoms;
  
 /**
  * Constructs a molecule and all atoms within it.  Atoms are arranged in
  * linked list that chains the "first" atom of the molecule to the "last".
  * @see Atom
  *
  * @param parent      species in which this molecule resides
  * @param n           number of atoms in this molecule
  */
  public Molecule(Species parent, int n) {
    parentSpecies = parent;
    Space.uEa1(r,0.0);
    Space.uEa1(p,0.0);
    zeroForce();
    neighbors = new Vector();
    setNAtoms(n);
  }
  
 /**
  * Constructs the molecules' atoms.  The number of atoms constructed
  * is equal to the current value of nAtoms.
  * Override this method if subclassing Molecule to define a molecule
  * comprising complex atoms (for example, nonspherical "atoms").
  */
/*  protected void makeAtoms() {
    atom = new Atom[nAtoms];
    for(int i=0; i<nAtoms; i++) {atom[i] = new Atom(this,i);}
  }
*/ 
 /**
  * Constructs atoms of molecule and links them, setting values of firstAtom and lastAtom.  
  * The number of atoms constructed is equal to the current value of nAtoms.
  * Does not handle linking to previous or next molecules' atoms
  */
  private final void makeAtoms() {
    Atom generator = parentSpecies.atomGenerator;
    firstAtom = generator.makeAtom(this,0);
    lastAtom = firstAtom;
    for(int i=1; i<nAtoms; i++) {
        lastAtom.setNextAtom(generator.makeAtom(this,i));
        lastAtom = lastAtom.getNextAtom();
    }
  }
    
 /**
  * @return the number of atoms in this molecule
  */
  public final int getNAtoms() {return nAtoms;}
  
 /**
  * Sets the number of atoms in the molecule, makes and links them, and
  * sets the molecular mass.  Does not handle linking of atoms to atoms
  * of other molecules.
  */
  private final void setNAtoms(int n) {
    nAtoms = n;
    makeAtoms();
//    orderAtoms();
    updateMass();
  }

  /**
   * Performs the list-linking of atoms in this molecule from the
   * firstAtom to the lastAtom.  Does not handle linking of atoms
   * to atoms of nextMolecule and previousMolecule
   * 
   * @see #setNextMolecule
   */
 /* private final void orderAtoms() {
    firstAtom = atom[0];
    lastAtom = atom[nAtoms-1];
    for(int i=1; i<nAtoms; i++) {atom[i-1].setNextAtom(atom[i]);}
  }    
  */
 /**
  * @return the next molecule following this one in the linked list of molecules,
  *         or null, if this is the lastMolecule in its phase
  */
  public final Molecule getNextMolecule() {return nextMolecule;}
  
 /**
  * Sets the argument to be this molecule's nextMolecule.  Also sets
  * the argument's previousMolecule to be this molecule, and does the
  * correcting linking of the respective molecules' first and last atoms.
  *
  * @param m   The molecule (possibly null) to be identified as this molecule's nextMolecule.
  */
  public final void setNextMolecule(Molecule m) {
    this.nextMolecule = m;
    if(m == null) {
        this.lastAtom.setNextAtom(null);
        return;
    }
    m.previousMolecule = this;
    this.lastAtom.setNextAtom(m.firstAtom);
  }
  
  public final void clearPreviousMolecule() {  //use setNextMolecule to set previousMolecule to non-null
    previousMolecule = null;
    firstAtom.clearPreviousAtom();
  }
 
 /**
  * @return the previous molecule before this one in the linked list of molecules,
  *         or null, if this is the firstMolecule in its phase
  */
  public final Molecule getPreviousMolecule() {return previousMolecule;}  
  
 /** Returns the squared distance of the center-of-mass of this
   * molecule relative to its value since the last call to updateNeighborList
   * Use is discouraged
   */
  public final double getSquareDisplacement() {return parentSpecies.parentPhase.space.r1Mr2_S(COM(),r0);}
  
 /**
   * Returns the intramolecular potential that governs interactions of all 
   * atoms of this molecule.
   *
   * @return  single-molecule potential for this molecule's atoms
   * @see Potential1
   */
  public final Potential1 getP1() {return parentSpecies.parentPhase.potential1[getSpeciesIndex()];}
 
 /**
  * Returns a vector of intermolecular potentials between this molecule's
  * species and all other species in the system.
  *
  * @return  vector of intermolecular potentials between this molecule and all other molecules in the system.
  * @see Potential2
  */
  public final Potential2[] getP2() {return parentSpecies.parentPhase.potential2[getSpeciesIndex()];}
 
 /**
  * Sets to zero this molecule's force vector
  */
  public final void zeroForce() {Space.uEa1(f, 0.0);}
  
 /**
  * Increments this molecule's force vector
  *
  * @param force the force vector being added to this molecule's force
  */
  public final void addForce(double[] force) { Space.uPEv1(f, force);}
 /**
  * Decrements this molecule's force vector
  *
  * @param force the force vector being subtracted from this molecule's force
  */
  public final void subtractForce(double[] force) {Space.uMEv1(f, force);}  
  
 /**
  * Each species has a unique integer index that is used to identify the correct
  * intra- and inter-molecular potentials for its molecules.
  *
  * @return  the index of this molecule's parent species
  */
 public final int getSpeciesIndex() {return parentSpecies.getSpeciesIndex();}
  
 /**
  * @return  the parent species of this molecule
  */
  public final Species getSpecies() {return parentSpecies;}
 
 /**
  * @return the phase in which this molecule resides
  */
  public final Phase getPhase() {return parentSpecies.parentPhase;}
 
 /**
  * Returns stored value of this molecule's mass; assumes that atom masses and/or number
  * of atoms is same as when updateMass was last invoked
  * @see #updateMass
  * @return the total mass of this molecule, in amu
  */
  public final double getMass() {return this.mass;}
  
 /**
  * Computes the total mass of this molecule by summing the masses
  * of its constituent atoms.  Also updates the center-of-mass fraction
  * of all its atoms.
  *
  * @see Atom#updateCOMFraction
  */
  public final void updateMass() {
    this.mass = 0.0;
    for(Atom a=firstAtom(); a!=terminationAtom(); a=a.getNextAtom()) {this.mass += ((AtomC)a).getMass();}
    for(Atom a=firstAtom(); a!=terminationAtom(); a=a.getNextAtom()) {((AtomC)a).updateCOMFraction();}
  }
  
 /**
  * Computes and returns the total (intra- and inter-molecular) contribution
  * of this molecule to the system's potential energy.
  *
  * @return  this molecule's intra- and inter-molecular potential energy, divided by kB, in Kelvin
  */
  public final double potentialEnergy() {
    double energy = 0.0;
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    if(nAtoms > 1) {
        for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {
            energy += a.intraPotentialEnergy();
        }
        energy *= 0.5;  //remove double-counting
    }
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {
        energy += a.interPotentialEnergy();
    }
    return energy;
  }
  
 /**
  * Computes and returns this molecule's total kinetic energy (including
  * contributions from all internal motions)
  *
  * @return  kinetic energy in (amu)(Angstrom)<sup>2</sup>(ps)<sup>-2</sup>
  */
  public double kineticEnergy() {
    double energy = 0.0;
    Atom nextMoleculeAtom = lastAtom.getNextAtom();  //so not computed each time through loop
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {energy += ((AtomC)a).kineticEnergy();}
    return energy;
  }
  
    /**
     * Displaces all atoms in this molecule by a vector dr.  Does not store original coordinates.
     *
     * @param dr  vector specifying change in position
     * @see #displace
     */
  public final void translate(double[] dr) {
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {((AtomC)a).translate(dr);}
  }
    /**
     * Displaces all atoms in this molecule the distance dr in one coordinate direction.
     *
     * @param i   index specifying coordinate direction of displacement
     * @param dr  displacement distance
     */
  public final void translate(int i, double dr) {
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {((AtomC)a).translate(i,dr);}
  }
  
    /**
     * Displaces all atoms in molecule by a vector dr, saving their original positions,
     * which can be recovered by call to replace.  Useful for moving molecules in 
     * Monte Carlo trials.
     *
     * @param dr   vector specifying change in position
     * @see #replace
     * @see #translate
     */
   public final void displace(double[] dr) {   //need to revise this
    Atom nextMoleculeAtom = terminationAtom();
    for(AtomC a=(AtomC)firstAtom; a!=nextMoleculeAtom; a=(AtomC)a.getNextAtom()) {a.displace(dr);}
   }
   
     /**
      * Puts molecule's atoms back in position held before last call to displace
      *
      * @see #displace
      */
   public final void replace() {    //need to revise this
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(AtomC a=(AtomC)firstAtom; a!=nextMoleculeAtom; a=(AtomC)a.getNextAtom()) {a.replace();}
   }

  /**
   * Rescales this particle's center of mass to new position corresponding to a 
   * scaling up or down of the size of the space it is in.  Used for moving molecules
   * when doing isobaric-simulation volume changes.  All atoms in this molecule
   * keep their relative distances and orientations.  Uses displace, so original coordinates may
   * be recovered by a call to replace
   *
   * @see replace
   */
   public void inflate(double scale) {
      Space.uEa1Tv1(dr,scale-1.0,COM());
      this.displace(dr);
   }
    
  /**
   * Computes and returns the center-of-mass vector of this molecule.  Value is not stored,
   * so it is computed afresh with each call
   *
   * @return center-of-mass coordinate vector of this molecule, in Angstroms
   */
  public final double[] COM() {
    if(nAtoms == 1) {return ((AtomC)firstAtom).r;}
    Space.uEa1(r,0.0);
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(AtomC a=(AtomC)firstAtom; a!=nextMoleculeAtom; a=(AtomC)a.getNextAtom()) {
        Space.uPEa1Tv1(r,a.getCOMFraction(),a.r);
    }
    return r;
  }
  
  public final void setCOM(double[] newCOM) {
    double[] oldCOM = new double[Space.D];
    Space.uEv1(oldCOM,COM());
    Space.uTEa1(oldCOM,-1.0);
    translate(oldCOM);  //zero center-of-mass
    translate(newCOM);  //move to new COM
  }
  
 /* public final void accelerate(int i, double dp) {
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {a.accelerate(i,dp);}
  }*/

/**
 * Used for neighbor listing; current structure is being deprecated and should not be used
  public final boolean needNeighborUpdate() {
    return (getSquareDisplacement() > parentSpecies.neighborUpdateSquareDisplacement);
  }
 */
  public final void addNeighbor(Molecule m) {neighbors.addElement(m);}
  public final void clearNeighborList() {neighbors.removeAllElements();}
  public final Enumeration getNeighborList() {return neighbors.elements();}
  public final boolean hasNeighbor(Molecule m) {return neighbors.contains(m);}
  public final void updateNeighborList() {
    clearNeighborList();
    Potential2[] p2 = parentSpecies.parentPhase.potential2[getSpeciesIndex()];
    for(Molecule m=previousMolecule; m!=null; m=m.getPreviousMolecule()) {
        if(!m.hasNeighbor(this)) {
            if(p2[m.getSpeciesIndex()].isNeighbor(this,m)) {m.addNeighbor(this);}
        }
    }
    for(Molecule m=nextMolecule; m!=null; m=m.getNextMolecule()) {
        if(p2[m.getSpeciesIndex()].isNeighbor(this,m)) {addNeighbor(m);}
    }
    Space.uEv1(r0,COM());
  }
  
  /**
   * Draws this molecule by calling its atoms' draw routines. 
   *
   * @param g         graphics object to which atom is drawn
   * @param origin    origin of drawing coordinates (pixels)
   * @param scale     factor determining size of drawn image relative to
   *                  nominal drawing size
   * @see Atom#draw
   */
  public void draw(Graphics g, int[] origin, double scale) {
    Atom terminator = terminationAtom();
    for(Atom a=firstAtom; a!=terminator; a=a.getNextAtom()) {a.draw(g, origin, scale);}
  }
  
  public final Atom firstAtom() {return firstAtom;}
  public final Atom lastAtom() {return lastAtom;}
  public final Atom terminationAtom() {return lastAtom.getNextAtom();}
  
}
