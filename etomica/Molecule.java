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
   *  Makes molecule with specified number of atoms, all of the same type.
   *  Number of atoms on molecule can be changed.
   */
  public Molecule(Species ps, Phase pp, AtomType type, int n) {  
    parentSpecies = ps;
    parentPhase = pp;
    coordinate = parentSpecies.parentSimulation.space.makeCoordinate();
    r = coordinate.position();
    p = coordinate.momentum();
    temp = parentPhase.parentSimulation.space.makeVector();
    nAtoms = n;
    
    firstAtom = new Atom(this,type,0);
    lastAtom = firstAtom;
    for(int i=1; i<nAtoms; i++) {
        lastAtom.setNextAtom(new Atom(this,type,i));
        lastAtom = lastAtom.nextAtom();
    }
  }
  
  /**
   *  Makes molecule with atoms of type specified in AtomType array.  Number of atoms is given by length of array. 
   *  Each atom may be of different type, or of same type but with different parameter values.
   *  Number of atoms cannot be changed.
   *  If all atoms are of same type with same parameters, can use the other constructor.
   */
  public Molecule(Species ps, Phase pp, AtomType[] type) {  
    parentSpecies = ps;
    parentPhase = pp;
    coordinate = parentSpecies.parentSimulation.space.makeCoordinate();
    r = coordinate.position();
    p = coordinate.momentum();
    nAtoms = type.length;
    temp = parentPhase.parentSimulation.space.makeVector();
    
    firstAtom = new Atom(this,type[0],0);
    lastAtom = firstAtom;
    for(int i=1; i<nAtoms; i++) {
        lastAtom.setNextAtom(new Atom(this,type[i],i));
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
   * Returns the intramolecular potential that governs interactions of all 
   * atoms of this molecule.
   *
   * @return  single-molecule potential for this molecule's atoms
   * @see Potential1
   */
  public final Potential1 getP1() {return parentSpecies.parentSimulation.potential1[getSpeciesIndex()];}
 
 /**
  * Returns a vector of intermolecular potentials between this molecule's
  * species and all other species in the system.
  *
  * @return  vector of intermolecular potentials between this molecule and all other molecules in the system.
  * @see Potential2
  */
  public final Potential2[] getP2() {return parentSpecies.parentSimulation.potential2[getSpeciesIndex()];}
   
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
  public final Phase getPhase() {return parentPhase;}
  
  public final double mass() {
    if(nAtoms==1) {return firstAtom.type.mass();}
    double mass = 0.0;
    for(Atom a=firstAtom(); a!=terminationAtom(); a=a.nextAtom()) {mass += a.type.mass();}
    return mass;
  }
  
  public final double rm() {
    if(nAtoms==1) {return firstAtom.type.rm();}
    else {return 1.0/mass();}
  }
  
 //  public void inflate(double scale) {coordinate.inflate(scale);}
/*      Space.uEa1Tv1(dr,scale-1.0,COM());
      this.displace(dr);
   }*/
    
  /**
   * Computes and returns the center-of-mass vector of this molecule.  Value is not stored,
   * so it is computed afresh with each call
   *
   * @return center-of-mass coordinate vector of this molecule, in Angstroms
   */
   public final Space.Vector COM() {
     updateR();
     return r;
   }
   public final void setCOM(Space.Vector u) {
     translateTo(u);
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
    for(Atom a=firstAtom; a!=terminator; a=a.nextAtom()) {a.draw(g, origin, scale);}
  }
  
  public final Atom firstAtom() {return firstAtom;}
  public final Atom lastAtom() {return lastAtom;}
  public final Atom terminationAtom() {return lastAtom.nextAtom();}
  
  Phase parentPhase;
  
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
  public final Space.Coordinate coordinate;   //might want private becuase coordinate must be evaluated from atom coordinate
  public final Space.Vector r;
  public final Space.Vector p;
//-------------------------------------------------------------------------

        protected final Space.Vector temp;
        public void updateR() {  //recomputes COM position from atom positions
            Space.AtomCoordinate c = firstAtom.coordinate;
            if(nAtoms==1) {r.E(c.position());}  //one atom in molecule
            else {  //multiatomic
                r.Ea1Tv1(c.atom().mass(),c.position());
                do {c=c.nextCoordinate(); r.PEa1Tv1(c.atom().mass(),c.position());} while (c.atom()!=lastAtom);
                r.DE(mass());
            }
        }
        public void updateP() {  //recomputes total momentum from atom momenta
            Space.AtomCoordinate c = firstAtom.coordinate;
            p.E(c.momentum());
            if(nAtoms==1) {return;}  //one atom in molecule
            do {c=c.nextCoordinate(); p.PE(c.momentum());} while (c.atom()!=lastAtom);
        }
        public void translateToward(Space.Vector u, double d) {temp.Ea1Tv1(d,u); translateBy(temp);}
        public void translateBy(Space.Vector u) {
            Space.AtomCoordinate c = firstAtom.coordinate;
            c.translateBy(u);
            if(nAtoms == 1) {return;}
            do {c=c.nextCoordinate(); c.translateBy(u);} while (c.atom()!=lastAtom);
        }
        public void translateTo(Space.Vector u) {
            updateR();  //update COM vector
            temp.E(u);  //temp = destination vector
            temp.ME(r);   //temp = destination - original = dr
            translateBy(temp);
        }
        public void displaceBy(Space.Vector u) {
            Space.AtomCoordinate c = firstAtom.coordinate;
            c.displaceBy(u);
            if(nAtoms == 1) {return;}
            do {c=c.nextCoordinate(); c.displaceBy(u);} while (c.atom()!=lastAtom);
        }
        public void displaceTo(Space.Vector u) {
            updateR();  //update COM vector
            temp.E(u);  //temp = destination vector
            temp.ME(r);   //temp = destination - original = dr
            displaceBy(temp);
        }
        public void displaceWithin(double d) {
            temp.setRandom(d);
            displaceBy(temp);
        }
            
        public void displaceToRandom(Phase phase) {displaceTo(phase.boundary().randomPosition());}
        public void translateToRandom(Phase phase) {translateTo(phase.boundary().randomPosition());}

        public void replace() {
            Space.AtomCoordinate c = firstAtom.coordinate;
            c.replace();
            if(nAtoms == 1) {return;}
            do {c=c.nextCoordinate(); c.replace();} while (c.atom()!=lastAtom);
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
        public void inflate(double s) {
            updateR();
            temp.Ea1Tv1(s-1.0,r);
            displaceBy(temp);   //displaceBy doesn't use temp
        }
        public Space.Vector position() {updateR(); return r;}
        public Space.Vector momentum() {updateP(); return p;}
        public double position(int i) {updateR(); return r.component(i);}
        public double momentum(int i) {updateP(); return p.component(i);}
 /**
  * Computes and returns this molecule's total kinetic energy (including
  * contributions from all internal motions)
  *
  * @return  kinetic energy in (amu)(Angstrom)<sup>2</sup>(ps)<sup>-2</sup>
  */
        public double kineticEnergy() {updateP(); return 0.5*p.squared()*rm();}
        
        public void randomizeMomentum(double temperature) {
            Space.AtomCoordinate c = firstAtom.coordinate;
            c.randomizeMomentum(temperature);
            if(nAtoms == 1) {return;}
            do {c=c.nextCoordinate(); c.randomizeMomentum(temperature);} while (c.atom()!=lastAtom);
        }

 /**
  * @return the previous molecule before this one in the linked list of molecules,
  *         or null, if this is the firstMolecule in its phase
  */
 public final Molecule previousMolecule() {return previousMolecule;}
 /**
  * @return the next molecule following this one in the linked list of molecules,
  *         or null, if this is the lastMolecule in its phase
  */
 public final Molecule nextMolecule() {return nextMolecule;}
  
  
 /**
  * Sets the argument to be this molecule's nextMolecule.  Also sets
  * the argument's previousMolecule to be this molecule, and does the
  * correcting linking of the respective molecules' first and last atoms.
  *
  * @param m   The molecule (possibly null) to be identified as this molecule's nextMolecule.
  */
  public final void setNextMolecule(Molecule m) {
    nextMolecule = m;
    if(m==null) {lastAtom.setNextAtom(null);}
    else {
        m.previousMolecule = this;
        lastAtom.setNextAtom(m.firstAtom);
    }
  }
  
  public final void clearPreviousMolecule() {  //use setNextMolecule to set previousMolecule to non-null
    previousMolecule = null;
    firstAtom.clearPreviousAtom();
  }
  

}
