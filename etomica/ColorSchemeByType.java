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
  
public final class ColorSchemeByType extends ColorScheme {
    
    public ColorSchemeByType() {}
        
 /**
    Return without changing atom's color
  */
    public final void setAtomColor(Atom a) {}
    
 /**
  * Initialize atom color to the color of its type
  */
    public final void initializeAtomColor(Atom a) {a.setColor(a.type.color());}

}
