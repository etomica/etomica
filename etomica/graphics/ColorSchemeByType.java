package etomica.graphics;
import etomica.*;
import java.awt.Color;

/**
 * Colors the atom according to the color given by its type field.  Instantiation
 * of this class leads to the placement of a ParameterColor object with each
 * existing and subsequent instance of AtomType.  This is used to reference the
 * color associated with each type.
 *
 * @author David Kofke
 */

public final class ColorSchemeByType extends ColorScheme implements Parameter.Source {
    
    public static final int colorIndex = AtomType.requestParameterIndex(
        new Parameter.Source() {
            public Parameter makeParameter() {return new ParameterColor();}
        });
        
    public ColorSchemeByType() {}
            
 /**
  * Initialize atom color to the color of its type
  */
    public final Color atomColor(Atom a) {
        return(((ParameterColor)a.type.parameter[colorIndex]).getColor());}
  
    //implementation of Parameter.Source interface
    public Parameter makeParameter() {return new ParameterColor();}
}
