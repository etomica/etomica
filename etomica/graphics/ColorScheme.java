package etomica.graphics;
import java.awt.Color;

import etomica.atom.AtomLeaf;

/**
 * Class that defines the algorithm used to determine atoms colors when drawn to DisplayPhase.
 * The atomColor method is called just before the atom is drawn to set the graphics color.
 *
 * @author David Kofke
 */
 
public abstract class ColorScheme implements java.io.Serializable {

    protected Color defaultColor;
    
    public ColorScheme() {
        this(DEFAULT_ATOM_COLOR);
    }
    public ColorScheme(Color color) {
        defaultColor = color;
    }
    
    public abstract Color getAtomColor(AtomLeaf a);
    
    public final void setDefaultColor(Color c) {defaultColor = c;}
    public final Color getDefaultColor() {return defaultColor;}

    public static Color DEFAULT_ATOM_COLOR = Color.red;
    
    /**
     * Colors all atoms with baseColor.
     */
    public static class Simple extends ColorScheme {
        public Simple() {super();}
        public Simple(java.awt.Color color) {super(color);}
        public Color getAtomColor(AtomLeaf a) {return defaultColor;}
    }//end of Simple
}//end of ColorScheme
