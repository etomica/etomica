package simulate;
import java.io.*;
import java.util.*;
import java.awt.Graphics;

public class Molecule implements Serializable {

  int speciesIndex;
  double mass;
  Vector neighbors;
  private final double[] r0 = new double[Space.D];
  private final double[] r = new double[Space.D];
  private final double[] dr = new double[Space.D];
  final double[] p = new double[Space.D];  //maybe delete p and f some day
  final double[] f = new double[Space.D];
  Species parentSpecies;
  Molecule nextMolecule, previousMolecule;
  Atom firstAtom, lastAtom;
  Atom[] atom;
  int nAtoms;
  
  public Molecule(Species s, int n) {
    parentSpecies = s;
    Space.uEa1(r,0.0);
    Space.uEa1(p,0.0);
    zeroForce();
    neighbors = new Vector();
    setNAtoms(n);
  }
  
  //Override for unconventional atoms
  protected void makeAtoms() {
    atom = new Atom[nAtoms];
    for(int i=0; i<nAtoms; i++) {atom[i] = new Atom(this,i);}
  }
  
  public final int getNAtoms() {return nAtoms;}
  private final void setNAtoms(int n) {
    nAtoms = n;
    makeAtoms();
    orderAtoms();
    updateMass();
  }
  
  private final void orderAtoms() {
    firstAtom = atom[0];
    lastAtom = atom[nAtoms-1];
    for(int i=1; i<nAtoms; i++) {atom[i-1].setNextAtom(atom[i]);}
  }    
  
  public final Molecule getNextMolecule() {return nextMolecule;}
  public final void setNextMolecule(Molecule m) {
    this.nextMolecule = m;
    if(m == null) {
        this.lastAtom.setNextAtom(null);
        return;
    }
    m.previousMolecule = this;
    this.lastAtom.setNextAtom(m.firstAtom);
  }
  public final Molecule getPreviousMolecule() {return previousMolecule;}  
  
  public final double getSquareDisplacement() {return parentSpecies.parentPhase.space.r1Mr2_S(COM(),r0);}
  
  public final Potential1 getP1() {return parentSpecies.parentPhase.potential1[getSpeciesIndex()];}
  public final Potential2[] getP2() {return parentSpecies.parentPhase.potential2[getSpeciesIndex()];}
    
  public final void zeroForce() {Space.uEa1(f, 0.0);}
  public final void addForce(double[] force) { Space.uPEv1(f, force);}
  public final void subtractForce(double[] force) {Space.uMEv1(f, force);}  
  
  public final int getSpeciesIndex() {return parentSpecies.getSpeciesIndex();}
  
  public final Species getSpecies() {return parentSpecies;}
  public final Phase getPhase() {return parentSpecies.parentPhase;}
  
  public final double getMass() {return this.mass;}
  public final void updateMass() {
    this.mass = 0.0;
    for(int i=0; i<nAtoms; i++) {this.mass += atom[i].getMass();}
    for(int i=0; i<nAtoms; i++) {atom[i].updateCOMFraction();}
  }
  
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
  
  public double kineticEnergy() {
    double energy = 0.0;
    Atom nextMoleculeAtom = lastAtom.getNextAtom();  //so not computed each time through loop
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {energy += a.kineticEnergy();}
    return energy;
  }
  
  public final void translate(double[] dr) {
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {a.translate(dr);}
  }
  public final void translate(int i, double dr) {
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {a.translate(i,dr);}
  }
  
   public final void displace(double[] dr) {
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {a.displace(dr);}
   }
   public final void replace() { 
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {a.replace();}
   }
   
   public final void inflate(double scale) {
      Space.uEa1Tv1(dr,scale-1.0,COM());
      this.displace(dr);
   }
    
  
  public final double[] COM() {
    if(nAtoms == 1) {return firstAtom.r;}
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

  public final boolean needNeighborUpdate() {
    return (getSquareDisplacement() > parentSpecies.neighborUpdateSquareDisplacement);
  }
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
  
  public void draw(Graphics g, int[] origin, double scale) {
    Atom nextMoleculeAtom = lastAtom.getNextAtom();
    for(Atom a=firstAtom; a!=nextMoleculeAtom; a=a.getNextAtom()) {a.draw(g, origin, scale);}
  }
  
}
