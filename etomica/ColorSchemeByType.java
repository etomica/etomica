package etomica;
import java.awt.Color;
/**
 * Colors the atom according to the color given by its type field.
 * @author David Kofke
 *
 */

public final class ColorSchemeByType extends ColorScheme {
    
    public ColorSchemeByType() {}
            
 /**
  * Initialize atom color to the color of its type
  */
    public final void colorAtom(Atom a) {a.setColor(a.type.color());}
  
}
