package etomica.graphics;
import java.awt.Color;

import etomica.Atom;

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
    public final Color atomColor(Atom a) {return null;}

}
