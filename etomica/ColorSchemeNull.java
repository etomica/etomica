package simulate;
import java.awt.Color;
/**
 *
 * @author David Kofke
 *
 */

 /**
  * Does nothing at any time to set atom's color, leaving color to be set type AtomType
  */
  
public final class ColorSchemeNull extends ColorScheme {
    
    public ColorSchemeNull() {}
        
 /**
    Return without changing atom's color
  */
    public final void setAtomColor(Atom a) {}
    public final void initializeAtomColor(Atom a) {}

}
