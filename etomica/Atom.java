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
public class Atom {

    /**
     * Constructs an atom with default values for mass, diameter, and color.
     * Default values are mass = 1.0 amu; diameter = 0.1 A; color = black.
     * Defaults for all coordinates and momenta are zero.
     *
     * @param parent       molecule in which atom resides
     * @param index        sequential index of atom as assigned by parent molecule
     */
    public Atom(Molecule parent, AtomType t, int index) {
        parentMolecule = parent;
        type = t;
        atomIndex = index;
        coordinate = parentMolecule.parentSpecies.parentSimulation.space.makeAtomCoordinate(this);
        setStationary(false);
        useTypeColor();
    }
                
    public void setIntegratorAgent(Integrator.Agent ia) {this.ia = ia;}
    
    public final Molecule getMolecule() {return parentMolecule;}
    
    public final int getSpeciesIndex() {return parentMolecule.getSpeciesIndex();}
    public final int getAtomIndex() {return atomIndex;}
    
    public final Color getColor() {return color;}
    public final void setColor(Color c) {this.color = c;}
    public final void useTypeColor() {this.color = type.color();}  //indicates that atom color is determined by its type
    
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

  public final double mass() {return type.mass();}
  public final double rm() {return type.rm();}
  
  public final Phase phase() {return parentMolecule.parentPhase;}

  public void draw(Graphics g, int[] origin, double scale) {type.draw(g, origin, scale, color, coordinate);}

    public Integrator.Agent ia;
    
    /**
     * Color of the atom when drawn on the screen
     * This color is set by the colorScheme object in the atom's species
     */
    Color color = Color.black;
    
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
    
    public final Space.AtomCoordinate coordinate;
    
    public final AtomType type;
    
}