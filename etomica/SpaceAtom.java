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
public abstract class SpaceAtom {

    /**
     * Constructs an atom with no initialization if parent is null; otherwise constructs atom with default atomIndex = 0.  
     * Expected use of such an Atom is for the construction of other Atoms via makeAtom method
     */
    public SpaceAtom(Molecule parent) {   //delete this when possible
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
    public SpaceAtom(Molecule parent, int index) {
        parentMolecule = parent;
        atomIndex = index;
        setStationary(false);
        setColor(Color.black);
    }
        
//    public abstract SpaceAtom makeAtom(Molecule m, int i);
    
    public void setIntegratorAgent(IntegratorAgent ia) {this.ia = ia;}
    
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
   * Sets atom following this one in linked list, and sets this to be that
   * atom's previous atom in list
   * 
   * @param atom  the next atom in the list
   * @see Molecule#orderAtoms
   */
  public void setNextAtom(SpaceAtom atom) {
    this.nextAtom = atom;
    if(atom != null) {atom.previousAtom = this;}
  }
  public void clearPreviousAtom() {previousAtom = null;}
  
  public final SpaceAtom nextMoleculeFirstAtom() {return parentMolecule.lastAtom.nextAtom;}  //first atom on next molecule
  public final SpaceAtom previousMoleculeLastAtom() {return parentMolecule.firstAtom.previousAtom;}  //first atom on next molecule

  public abstract void draw(Graphics g, int[] origin, double scale);

//    public IntegratorHard integrator = new IntegratorHard();
//    public IntegratorAgent ia = integrator.makeAgent(this);
    public IntegratorAgent ia;
    
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
     * next atom in linked list of atoms.  
     * Because it may be shadowed in subclasses, it should be set only via call to setNextAtom.
     * @see #setNextAtom
     */
    public SpaceAtom nextAtom = null;
    
    /**
     * previous atom in linked list of atoms.  
     * Becasue it may be shadowed in subclasses, it should be set only via call to setPreviousAtom or clearPreviousAtom.
     *
     * @see #setNextAtom
     */
    public SpaceAtom previousAtom = null;
        
    /**
     * Flag indicating whether atom is stationary or mobile.
     * Default is false (atom is mobile)
     */
    private boolean stationary;
   
}