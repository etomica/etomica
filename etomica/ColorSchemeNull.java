package simulate;
import java.awt.Color;

/**
 * Does nothing at any time to set atom's color, leaving color to be set to default value.
 * @author David Kofke
 *
 */

public final class ColorSchemeNull extends ColorScheme {
    
    public ColorSchemeNull() {}
        
 /**
  * Return without changing atom's color.
  */
    public final void colorAtom(Atom a) {}
 /**
  * Override to remove loop over atoms.
  */
    public final void colorAllAtoms(Atom a) {}

}
