package etomica;
import java.util.EventObject;

/**
 * Event that conveys some activity with respect to a phase or the things it contains.
 *
 * @see PhaseListener
 * @see DisplayPhaseListener
 */
public class PhaseEvent extends java.util.EventObject {
    
    protected Phase phase;
    protected Molecule molecule;
    protected Atom atom;
    protected Space.Vector point;
    
    public PhaseEvent(Object source) {
        super(source);
    }
    public static int POINT_SELECTED = 0;
    public static int ATOM_SELECTED = 10;
    public static int ATOM_RELEASED = 11;
    public static int MOLECULE_SELECTED = 20;
    public static int MOLECULE_RELEASED = 21;
    
    public void setPhase(Phase p) {phase = p;}
    public Phase getPhase() {return phase;}
    
    public void setPoint(Space.Vector p) {point = p;}
    public Space.Vector getPoint() {return point;}
    
    public void setAtom(Atom a) {atom = a;}
    public Atom getAtom() {return atom;}
    
    public void setMolecule(Molecule m) {molecule = m;}
    public Molecule getMolecule() {return molecule;}
}
    