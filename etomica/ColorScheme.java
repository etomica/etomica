package etomica;
import java.awt.*;

/**
 * Class that defines the algorithm used to determine atoms colors when drawn to DisplayPhase.
 * The color of each atom when drawn is given by its color field.  The colorAllAtoms method of
 * this class is called just before the atoms are drawn; this method sets the color field of
 * each atom according to the coloring algorithm.
 *
 * @author David Kofke
 */
 
public abstract class ColorScheme implements java.io.Serializable {

    public static String getVersion() {return "01.07.13";}
 
    private Phase phase;
    protected Color baseColor;
    protected AtomIterator atomIterator;
    protected AtomAction action = new Action();
    
    public ColorScheme() {
        this(Default.ATOM_COLOR);
    }
    public ColorScheme(Color color) {
        baseColor = color;
    }
    
    public void setPhase(Phase p) {
        phase = p;
        setIterator(p.iteratorFactory().makeAtomIterator());
    }
    
    public Phase getPhase() {return phase;}
    
    public void setIterator(AtomIterator iterator) {
        atomIterator = iterator;
    }
    public AtomIterator getIterator() {return atomIterator;}
        
    public AtomAction action() {return action;}
        
    public abstract void colorAtom(Atom a);
    
    /**
     * Colors all atoms using the colorAtom method.
     * Called repeatedly with each redraw of phase, so if color scheme never
     * causes atom color to change, override this method so that it performs no action.
     * Then atoms are not re-colored redundantly every time drawn.
     */
    public void colorAllAtoms() {
        if(atomIterator != null) atomIterator.allAtoms(action);
    }
    
    /**
     * Colors all atoms using the colorAtom method.
     * Called only when color scheme is added to phase or when the iterator is changed.
     * In the ColorScheme base class, this method performs exactly the same action as the
     * colorAllAtoms method, but in subclasses colorAllAtoms may be overridden to eliminate
     * its action; then atom color is set only by call to this method.
     */
    public final void initialColorAllAtoms() {
        if(atomIterator != null) atomIterator.allAtoms(action);
    }
    
    public final void setBaseColor(Color c) {baseColor = c;}
    public final Color getBaseColor() {return baseColor;}

    public class Action extends AtomAction {
        
        public void actionPerformed(Atom a) {
            colorAtom(a);
        }
    }
    
    /**
     * Colors all atoms with baseColor.
     */
    public static class Simple extends ColorScheme {
        public Simple() {super();}
        public Simple(java.awt.Color color) {super(color);}
        public void colorAtom(Atom a) {a.setColor(baseColor);}
    }
}
