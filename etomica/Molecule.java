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
  
 /**
  * Center-of-mass (COM) coordinate
  */
  public final PhaseSpace.MoleculeCoordinate coordinate;   //might want private becuase coordinate must be evaluated from atom coordinate
  
  public Molecule(Species parent, int n) {
    parentSpecies = parent;
    coordinate = parentSpecies.parentPhaseSpace.makeMoleculeCoordinate(this);
    nAtoms = n;
    makeAtoms();
  }
  
 /**
  * Constructs atoms of molecule and links them, setting values of firstAtom and lastAtom.  
  * The number of atoms constructed is equal to the current value of nAtoms.
  * Does not handle linking to previous or next molecules' atoms
  * Previously existing atoms are discarded
  */
  private final void makeAtoms() {
    firstAtom = new Atom(this,0);
    lastAtom = firstAtom;
    for(int i=1; i<nAtoms; i++) {
        lastAtom.setNextAtom(new Atom(this,i));
        lastAtom = lastAtom.nextAtom();
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
  * Any existing atoms are discarded
  */
/*  private final void setNAtoms(int n) {
    nAtoms = n;
    makeAtoms();
  } */

 /**
  * @return the next molecule following this one in the linked list of molecules,
  *         or null, if this is the lastMolecule in its phase
  */
  public final Molecule nextMolecule() {return coordinate.nextMolecule();}
  
 /**
  * @return the previous molecule before this one in the linked list of molecules,
  *         or null, if this is the firstMolecule in its phase
  */
  public final Molecule previousMolecule() {return coordinate.previousMolecule();}  
  
 /**
  * Sets the argument to be this molecule's nextMolecule.  Also sets
  * the argument's previousMolecule to be this molecule, and does the
  * correcting linking of the respective molecules' first and last atoms.
  *
  * @param m   The molecule (possibly null) to be identified as this molecule's nextMolecule.
  */
  public final void setNextMolecule(Molecule m) {
    if(m==null) {
        coordinate.setNextCoordinate(null);
        lastAtom.setNextAtom(null);
    }
    else {
        coordinate.setNextCoordinate(m.coordinate);
        lastAtom.setNextAtom(m.firstAtom);
    }
  }
  
  public final void clearPreviousMolecule() {  //use setNextMolecule to set previousMolecule to non-null
    coordinate.clearPreviousCoordinate();
    firstAtom.clearPreviousAtom();
  }
  
 /**
   * Returns the intramolecular potential that governs interactions of all 
   * atoms of this molecule.
   *
   * @return  single-molecule potential for this molecule's atoms
   * @see Potential1
   */
  public final Potential1 getP1() {return parentSpecies.parentPhaseSpace.potential1[getSpeciesIndex()];}
 
 /**
  * Returns a vector of intermolecular potentials between this molecule's
  * species and all other species in the system.
  *
  * @return  vector of intermolecular potentials between this molecule and all other molecules in the system.
  * @see Potential2
  */
  public final Potential2[] getP2() {return parentSpecies.parentPhaseSpace.potential2[getSpeciesIndex()];}
   
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
  public final PhaseSpace getPhaseSpace() {return parentSpecies.parentPhaseSpace;}
   
/* /**
  * Computes the total mass of this molecule by summing the masses
  * of its constituent atoms.  Also updates the center-of-mass fraction
  * of all its atoms.
  *
  * @see Atom#updateCOMFraction
  */
  /*
  public final void updateMass() {
    this.mass = 0.0;
    for(Atom a=firstAtom(); a!=terminationAtom(); a=a.getNextAtom()) {this.mass += ((AtomC)a).getMass();}
    for(Atom a=firstAtom(); a!=terminationAtom(); a=a.getNextAtom()) {((AtomC)a).updateCOMFraction();}
  }*/
  
  
 /**
  * Computes and returns this molecule's total kinetic energy (including
  * contributions from all internal motions)
  *
  * @return  kinetic energy in (amu)(Angstrom)<sup>2</sup>(ps)<sup>-2</sup>
  */
  public final double kineticEnergy() {return coordinate.kineticEnergy();}
  
    /**
     * Displaces all atoms in this molecule by a vector dr.  Does not store original coordinates.
     *
     * @param dr  vector specifying change in position
     * @see #displace
     */
  public final void translateTo(PhaseSpace.Vector v) {coordinate.translateTo(v);}
  public final void translateBy(PhaseSpace.Vector v) {coordinate.translateBy(v);}

    /**
     * Displaces all atoms in molecule by a vector dr, saving their original positions,
     * which can be recovered by call to replace.  Useful for moving molecules in 
     * Monte Carlo trials.
     *
     * @param dr   vector specifying change in position
     * @see #replace
     * @see #translate
     */
   public final void displaceTo(PhaseSpace.Vector v) {coordinate.displaceTo(v);}
     /**
      * Puts molecule's atoms back in position held before last call to displace
      *
      * @see #displace
      */
   public final void replace() {coordinate.replace();}

  /**
   * Rescales this particle's center of mass to new position corresponding to a 
   * scaling up or down of the size of the space it is in.  Used for moving molecules
   * when doing isobaric-simulation volume changes.  All atoms in this molecule
   * keep their relative distances and orientations.  Uses displace, so original coordinates may
   * be recovered by a call to replace
   *
   * @see replace
   */
   public void inflate(double scale) {coordinate.inflate(scale);}
/*      Space.uEa1Tv1(dr,scale-1.0,COM());
      this.displace(dr);
   }*/
    
  /**
   * Computes and returns the center-of-mass vector of this molecule.  Value is not stored,
   * so it is computed afresh with each call
   *
   * @return center-of-mass coordinate vector of this molecule, in Angstroms
   */
   public final PhaseSpace.Vector COM() {
     coordinate.update();
     return coordinate.position();
   }
/*    if(nAtoms == 1) {return firstAtom.r;}
    Space.uEa1(r,0.0);
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {
        Space.uPEa1Tv1(r,a.getCOMFraction(),a.r);
    }
    return r;
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
/*  public final void addNeighbor(Molecule m) {neighbors.addElement(m);}
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
  }*/
  
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
    for(Atom a=firstAtom; a!=terminator; a=a.nextAtom()) {a.draw(g, origin, scale);}
  }
  
  public final Atom firstAtom() {return firstAtom;}
  public final Atom lastAtom() {return lastAtom;}
  public final Atom terminationAtom() {return lastAtom.nextAtom();}
  
}
