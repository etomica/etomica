package simulate;
import java.awt.*;

/**
 *  Each instance of the class Atom holds the position and velocity of one
 *  physical atom; all simulation kinetics and dynamics are performed by 
 *  operating on these values.
 *  
 *  @author David Kofke
 *  @author C. Daniel Barnes
 *  @see Molecule
 */
public abstract class Atom {

    /**
     * Color of the atom when drawn on the screen
     */
    Color color = Color.black;
    
    /**
     * Ratio of atom mass to mass of parent molecule.
     * Used to locate the center-of-mass of the molecule
     * @see Molecule#COM
     */
    double COMFraction = 1.0;

    /**
     * Instance of molecule in which this atom resides.
     * Assigned in Atom constructor.
     * @see Molecule#makeAtoms
     */
    Molecule parentMolecule = null;
    
    /**
     * Identifier of atom within molecule.
     * Assigned by parent molecule when invoking Atom constructor.
     * @see Molecule#makeAtoms
     */
    int atomIndex;
    
    
    /**
     * next atom in linked list of atoms
     * @see #setNextAtom
     */
    Atom nextAtom = null;
    
    /**
     * previous atom in linked list of atoms
     *
     * @see #setNextAtom
     */
    Atom previousAtom = null;
        
    /**
     * Flag indicating whether atom is stationary or mobile.
     * Default is false (atom is mobile)
     */
    private boolean stationary;
    
    /**
     * Constructs an atom with no initialization if parent is null; otherwise constructs atom with default atomIndex = 0.  
     * Expected use of such an Atom is for the construction of other Atoms via makeAtom method
     */
    public Atom(Molecule parent) {
        this(parent,0);
    }
    
    /**
     * Constructs an atom with default values for mass, diameter, and color.
     * Default values are mass = 1.0 amu; diameter = 0.1 A; color = black.
     * Defaults for all coordinates and momenta are zero.
     *
     * @param parent       molecule in which atom resides
     * @param index        sequential index of atom as assigned by parent molecule
     */
    public Atom(Molecule parent, int index) {
        parentMolecule = parent;
        atomIndex = index;
        setStationary(false);
        setColor(Color.black);
    }
        
    public abstract Atom makeAtom(Molecule m, int i);
    
    public final Molecule getMolecule() {return parentMolecule;}
    
    public final int getSpeciesIndex() {return parentMolecule.getSpeciesIndex();}
    public final int getAtomIndex() {return atomIndex;}
    
    public final Color getColor() {return color;}
    public final void setColor(Color c) {this.color = c;}
    
    public void setStationary(boolean b) {
        stationary = b;
    }
    public final boolean isStationary() {return stationary;}

    /**
     * Computes and returns the potential energy of the atom due to its interactions
     * with all other atoms (intra and intermolecular) in the phase
     *
     * @return (potential energy)/kB in Kelvins
     */
    public final double potentialEnergy() {
        return intraPotentialEnergy() + interPotentialEnergy();
    }
    
  public abstract double interPotentialEnergy();
  public abstract double intraPotentialEnergy();
    
  public final Atom getNextAtom() {return nextAtom;}
  /**
   * Sets atom following this one in linked list, and sets this to be that
   * atom's previous atom in list
   * 
   * @param atom  the next atom in the list
   * @see Molecule#orderAtoms
   */
  public void setNextAtom(Atom atom) {
    this.nextAtom = atom;
    if(atom != null) {atom.previousAtom = this;}
  }
  
  public final Atom nextMoleculeFirstAtom() {return parentMolecule.lastAtom.getNextAtom();}  //first atom on next molecule

  public final void clearPreviousAtom() {previousAtom = null;}
  
  public final Atom getPreviousAtom() {return previousAtom;}
  
  public abstract void draw(Graphics g, int[] origin, double scale);
}