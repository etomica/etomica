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
        coordinate = parentMolecule.parentSpecies.parentPhaseSpace.makeAtomCoordinate(this);
        setStationary(false);
        setColor(Color.black);
    }
                
    public void setIntegratorAgent(IntegratorAgent ia) {this.ia = ia;}
    
    public final Molecule getMolecule() {return parentMolecule;}
    
    public final int getSpeciesIndex() {return parentMolecule.getSpeciesIndex();}
    public final int getAtomIndex() {return atomIndex;}
    
    public final Color getColor() {return color;}
    public final void setColor(Color c) {this.color = c;}
    
    public void setStationary(boolean b) {stationary = b;}
    public final boolean isStationary() {return stationary;}

  /**
   * Sets atom following this one in linked list, and sets this to be that
   * atom's previous atom in list
   * 
   * @param atom  the next atom in the list
   */
  public final void setNextAtom(Atom atom) {
    if(atom==null) {coordinate.setNextCoordinate(null);}
    else {coordinate.setNextCoordinate(atom.coordinate);}
  }
  public final void clearPreviousAtom() {coordinate.clearPreviousCoordinate();}
  public final Atom nextAtom() {return coordinate.nextAtom();}
  public final Atom previousAtom() {return coordinate.previousAtom();}
  
  public final Atom nextMoleculeFirstAtom() {return parentMolecule.lastAtom.nextAtom();}  //first atom on next molecule
  public final Atom previousMoleculeLastAtom() {return parentMolecule.firstAtom.previousAtom();}  //first atom on next molecule

  public abstract void draw(Graphics g, int[] origin, double scale);

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
     * Flag indicating whether atom is stationary or mobile.
     * Default is false (atom is mobile)
     */
    private boolean stationary;
   
    /**
     * Instance of molecule in which this atom resides.
     * Assigned in Atom constructor.
     * @see Molecule#makeAtoms
     */
    final Molecule parentMolecule;
    
    /**
     * Identifier of atom within molecule.
     * Assigned by parent molecule when invoking Atom constructor.
     * @see Molecule#makeAtoms
     */
    final int atomIndex;
    
    public final PhaseSpace.AtomCoordinate coordinate;
}