package etomica;
import java.io.*;
import java.awt.Graphics;

//Java2 imports
//import java.util.Iterator;

import etomica.utility.Iterator;


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

public class Molecule implements Space.Occupant, Serializable {

  /**
   *  Makes molecule with specified number of atoms, all of the same type.
   *  Number of atoms on molecule can be changed, but need to modify code to take care of atomIterator if this is done.
   */
  public Molecule(Species ps, Phase pp, AtomType type, int n) {  
    parentSpecies = ps;
    setParentPhase(pp);
    coordinate = parentSimulation().space().makeCoordinate(this);
    r = coordinate.position();
    p = coordinate.momentum();
    temp = parentSimulation().space().makeVector();
    atomCount = n;
//    atomIterator = (atomCount > 1) ? new AtomIterator() : new MonoAtomIterator();
    if(atomCount > 1) atomIterator = new AtomIterator();
    else atomIterator = new MonoAtomIterator();
    
    
    firstAtom = new Atom(this,type,0);
    lastAtom = firstAtom;
    for(int i=1; i<atomCount; i++) {
        lastAtom.setNextAtom(new Atom(this,type,i));
        lastAtom = lastAtom.nextAtom();
    }
    parentSpecies.moleculeConfiguration.initializeCoordinates(this);
  }
  
  /**
   *  Makes molecule with atoms of type specified in AtomType array.  Number of atoms is given by length of array. 
   *  Each atom may be of different type, or of same type but with different parameter values.
   *  Number of atoms cannot be changed.
   *  If all atoms are of same type with same parameters, can use the other constructor.
   *  If mono-atomic and number of atoms not expected to change, may be advantageous to use this constructor
   */
  public Molecule(Species ps, Phase pp, AtomType[] type) {  
    parentSpecies = ps;
    setParentPhase(pp);
    atomCount = type.length;
    firstAtom = new Atom(this,type[0],0);
    lastAtom = firstAtom;
    for(int i=1; i<atomCount; i++) {
        lastAtom.setNextAtom(new Atom(this,type[i],i));
        lastAtom = lastAtom.nextAtom();
    }
    if(atomCount > 1) { //Multi-atomic
        atomIterator = new AtomIterator();
        coordinate = parentSimulation().space().makeCoordinate(this);
    }
    else {  //Mono-atomic
        atomIterator = new MonoAtomIterator();
        coordinate = firstAtom.coordinate();
    }
    r = coordinate.position();
    p = coordinate.momentum();
    temp = parentSimulation().space().makeVector();
    parentSpecies.moleculeConfiguration.initializeCoordinates(this);
    
  }
  
  public final Simulation parentSimulation() {return parentSpecies.parentSimulation();}
  public final Phase parentPhase() {return parentPhase;}
  public final Space.Coordinate coordinate() {
    updateR();
    updateP();
    return coordinate;
  }

 /**
  * @return the number of atoms in this molecule
  */
  public final int atomCount() {return atomCount;}
  
  /**
   * @return the ith atom of the molecule, with first at i=0
   */
   public Atom getAtom(int i) {
     if(i >= atomCount) throw new IllegalArgumentException("Request atom index greater than number of atoms in molecule");
     Atom a = firstAtom;
     for(int j=0; j<i; j++) {
        a = a.nextAtom();
     }
     return a;
   }
   
   /**
    * Integer index that can be set and accessed for any purpose, for example
    * to associate the molecule with a row/column of a matrix.
    */
    public void setIndex(int i) {index = i;}
   /**
    * Integer index that can be set and accessed for any purpose, for example
    * to associate the molecule with a row/column of a matrix.
    */
    public int getIndex() {return index;}

   
 /**
   * Returns the intramolecular potential that governs interactions of all 
   * atoms of this molecule.
   *
   * @return  single-molecule potential for this molecule's atoms
   * @see Potential1
   */
//  public final Potential1 getP1() {return potential1[speciesIndex()];}
 
 /**
  * Returns a vector of intermolecular potentials between this molecule's
  * species and all other species in the system.
  *
  * @return  vector of intermolecular potentials between this molecule and all other molecules in the system.
  * @see Potential2
  */
//  public final Potential2[] getP2() {return Simulation.potential2[speciesIndex()];}
   
 /**
  * Each species has a unique integer index that is used to identify the correct
  * intra- and inter-molecular potentials for its molecules.
  *
  * @return  the index of this molecule's parent species
  */
 public final int speciesIndex() {return parentSpecies.getSpeciesIndex();}
  
 /**
  * @return  the parent species of this molecule
  */
  public final Species parentSpecies() {return parentSpecies;}
 
 /**
  * @return the phase in which this molecule resides
  */
//  public final Phase getPhase() {return parentPhase;}
  
  public final double mass() {
    if(atomCount==1) {return firstAtom.type.mass();}
    double mass = 0.0;
    atomIterator.reset();
    while(atomIterator.hasNext()) {mass += atomIterator.next().type.mass();}
    return mass;
  }
  
  public final double rm() {
    if(atomCount==1) {return firstAtom.type.rm();}
    else {return 1.0/mass();}
  }
  
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
    atomIterator.reset();
    while(atomIterator.hasNext()) {atomIterator.next().draw(g, origin, scale);}
  }
  
  public final Atom firstAtom() {return firstAtom;}
  public final Atom lastAtom() {return lastAtom;}
  
  public final void setParentPhase(Phase p) {
    parentPhase = p;
    setContainer(p);
  }
  public final Molecule.Container container() {return container;}
  public final void setContainer(Molecule.Container mc) {container = mc;}
  
  private Phase parentPhase;
  
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
  int atomCount;
  
  private int index;
  
 /**
  * Center-of-mass (COM) coordinate
  */
  public final Space.Coordinate coordinate;   //might want private becuase coordinate must be evaluated from atom coordinate
  public final Space.Vector r;
  public final Space.Vector p;
  
  public  etomica.AtomIterator atomIterator;
  
  public Molecule.Container container;
  
    protected final Space.Vector temp;
            
    public void updateR() {  //recomputes COM position from atom positions
        if(atomCount==1) {r.E(firstAtom.coordinate().position());}  //one atom in molecule
        else {  //multiatomic
            r.E(0.0);
            atomIterator.reset();
            while(atomIterator.hasNext()) {
                Atom a = atomIterator.next();
                r.PEa1Tv1(a.mass(),a.coordinate().position());
            }
            r.DE(mass());
        }
    }
    public void updateP() {  //recomputes total momentum from atom momenta
        p.E(0.0);
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            p.PE(atomIterator.next().momentum());
        }
    }
    public void translateToward(Space.Vector u, double d) {temp.Ea1Tv1(d,u); translateBy(temp);}
    public void translateBy(Space.Vector u) {
        atomIterator.reset();
        while(atomIterator.hasNext()) {atomIterator.next().translateBy(u);}
    }
    public void accelerateBy(Space.Vector u) {
        atomIterator.reset();
        while(atomIterator.hasNext()) {atomIterator.next().accelerateBy(u);}
    }
    public void accelerateTo(Space.Vector u) {
        updateP();  //update COM vector
        temp.E(u);  //temp = destination vector
        temp.ME(p);   //temp = destination - original = dr
        accelerateBy(temp);
    }
    public void translateTo(Space.Vector u) {
        updateR();  //update COM vector
        temp.E(u);  //temp = destination vector
        temp.ME(r);   //temp = destination - original = dr
        translateBy(temp);
    }
    public void displaceBy(Space.Vector u) {
        atomIterator.reset();
        while(atomIterator.hasNext()) {atomIterator.next().displaceBy(u);}
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
//        public void translateToRandom(Phase phase) {translateTo(phase.boundary().randomPosition());}

    public void replace() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {atomIterator.next().replace();}
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
//        public void inflate(double s) {
//            updateR();
//            temp.Ea1Tv1(s-1.0,r);
//            displaceBy(temp);   //displaceBy doesn't use temp
//        }
        public Space.Vector position() {updateR(); return r;}
        public Space.Vector momentum() {updateP(); return p;}
        public double position(int i) {updateR(); return r.component(i);}
        public double momentum(int i) {updateP(); return p.component(i);}
 /**
  * Computes and returns this molecule's total kinetic energy (including
  * contributions from all internal motions)
  *
  * @return  kinetic energy in (Daltons)(Angstrom)<sup>2</sup>(ps)<sup>-2</sup>
  */
        public double kineticEnergy() {updateP(); return 0.5*p.squared()*rm();}
        
        public void randomizeMomentum(double temperature) {
            atomIterator.reset();
            while(atomIterator.hasNext()) {
                atomIterator.next().randomizeMomentum(temperature);
            }
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
  
  public final class AtomIterator extends etomica.AtomIterator {
        private Atom atom, nextAtom;
        private boolean hasNext;
        public AtomIterator() {reset();}
        public boolean hasNext() {return hasNext;}
        public Atom reset() {
            atom = firstAtom;
            hasNext = true;
            return atom;
        }
        public Atom reset(Atom a) {return reset();}
        public Atom reset(Atom a1, Atom a2) {return reset();}
        public Atom next() {
            nextAtom = atom;
            if(atom == lastAtom) {hasNext = false;}
            else {atom = atom.nextAtom();}
            return nextAtom;
        }
        public void allAtoms(AtomAction act) {
            Atom term = lastAtom.nextAtom();
            for(Atom a=firstAtom; a!=term; a=a.nextAtom()) {act.actionPerformed(a);}
        }
    } //end of AtomIterator
  public final class MonoAtomIterator extends etomica.AtomIterator {
        private boolean hasNext;
        public MonoAtomIterator() {reset();}
        public boolean hasNext() {return hasNext;}
        public Atom reset() {hasNext = true; return firstAtom;}
        public Atom reset(Atom a) {return reset();}
        public Atom reset(Atom a1, Atom a2) {return reset();}
        public Atom next() {
            hasNext = false;
            return firstAtom;
        }
        public void allAtoms(AtomAction act) {act.actionPerformed(firstAtom);}
    } //end of MonoAtomIterator

  /**
   * Molecule.Container
   * Interface for molecule reservoirs and phases
   */
  public interface Container {
 /** addMolecule should first remove molecule from original container, then
  *  do what is necessary to add molecule to new container
  */
      public void addMolecule(Molecule m);
  /**
   * removeMolecule should be called only by the addMolecule method of another container
   */
      public void removeMolecule(Molecule m);
  } //end of interface Container
  
/**
 * General class for assignment of coordinates to all atoms
 * Places atoms of each molecule in species in same position with respect to the molecule's center of mass
 */

    public static abstract class Configuration implements java.io.Serializable {
      
    protected Species parentSpecies;  //some subclasses may want to take an action on setting species, so don't make public
    protected final double[] dim;
      
//    public Configuration() {super();}
    public Configuration(Species parent) {
        parentSpecies = parent;
        dim = new double[parent.parentSimulation().space().D()];
    }
      
    public void initializeCoordinates() {
        if(parentSpecies == null) {return;}
        Iterator e = parentSpecies.agents.values().iterator();
        while(e.hasNext()) {
            Species.Agent agent = (Species.Agent)e.next();
            for(Molecule m=agent.firstMolecule(); m!=agent.terminationMolecule(); m=m.nextMolecule()) {
                initializeCoordinates(m);
            }
        }
        computeDimensions();
    }
            
            
    public void initializeCoordinates(Phase p) {
        if(parentSpecies == null) {return;}
        Species.Agent agent = parentSpecies.getAgent(p);
        for(Molecule m=agent.firstMolecule(); m!=agent.terminationMolecule(); m=m.nextMolecule()) {
            initializeCoordinates(m);
        }
        p.iteratorFactory().reset();
        computeDimensions();
    }
      
    public void setParentSpecies(Species s) {parentSpecies = s;}
    public Species parentSpecies() {return parentSpecies;}
      
    public double[] moleculeDimensions() {return dim;}
      
    public abstract void initializeCoordinates(Molecule m);
    protected abstract void computeDimensions();

       //Simple linear configuration of molecules
            public static class Linear extends Configuration {
                
                private double bondLength = 0.5*Default.ATOM_SIZE;
                private Space.Vector orientation;
                private double[] angle;
            //    private double theta = 45.;
                
//                public Linear(){
//                    orientation = Simulation.space.makeVector();
//                    setAngle(0,45.);
//                }
                public Linear(Species parent) {
//                    this();
                    super(parent);
                    orientation = parent.parentSimulation().space().makeVector();
                    angle = new double[parent.parentSimulation().space().D()];
                    setAngle(0,45.);
//                    parentSpecies = parent;
                }
              
                public void setBondLength(double b) {
                    bondLength = b;
                    computeDimensions();
                }
                public double getBondLength() {return bondLength;}
                
                public void setAngle(int i, double t) {
                    angle[i] = Math.PI*t/180.;
                    switch(parentSpecies.parentSimulation().space().D()) {
                        case 1:
                            return;
                        case 2:
                            setOrientation(new Space2D.Vector(Math.cos(angle[0]),Math.sin(angle[0])));
                            return;
            //            case 3:
            //                setOrientation(new Space3D.Vector(Math.sin(angle[1])*Math.cos(angle[0]),
            //                                                  Math.sin(angle[1])*Math.sin(angle[0]),
            //                                                  Math.cos(angle[1])));
            //                return;
                    }
                }
                public double getAngle(int i) {return angle[i];}
                public void setOrientation(Space.Vector e) {orientation.E(e);}
              
            /**
            * Sets all atoms coordinates to lie on a straight line along the x-axis, with the
            * center of mass unchanged from the value before method was called
            */
                public void initializeCoordinates(Molecule m) {
                    Space.Vector OldCOM = parentSpecies.parentSimulation().space().makeVector();
                    OldCOM.E(m.COM());
                    double xNext = 0.0;
                    m.atomIterator.reset();
                    while(m.atomIterator.hasNext()) {
                        Atom a = m.atomIterator.next();
                        a.translateTo(OldCOM);  //put all atoms at same point
                        a.translateBy(xNext,orientation);  //move xNext distance in direction orientation
                        xNext += bondLength;
                    }
                    m.translateTo(OldCOM);  //shift molecule to original COM
                }
                
                protected void computeDimensions() {
            /*        if(parentSpecies()==null) return;
                    dim[1] = 0.0;
                    Molecule m = parentSpecies().getMolecule();  //a typical molecule
                    initializeCoordinates(m);
                    for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.nextAtom()) {
                        dim[1] = Math.max(dim[1], ((AtomType.Disk)a.type).diameter());  //width is that of largest atom
                    }
                    dim[0] = 0.5*(((AtomType.Disk)m.firstAtom().type).diameter() + ((AtomType.Disk)m.lastAtom().type).diameter()) + (m.atomCount-1) * bondLength;
            */
                }

            public void setParentSpecies(Species s) {
                parentSpecies = s;
            //    setAngle(theta);   //fix orientation to x-axis of 2-D space
            //    computeDimensions();
            }
        } //end of Molecule.Configuration.Linear
    }//end of Molecule.Configuration
      
}
