package etomica;
import java.awt.Color;
/**
 * Sets color of a wall based on the value of its temperature field.
 * @author David Kofke
 *
 */
 
public class ColorSchemeWallTemperature extends ColorSchemeTemperature {
    
    public ColorSchemeWallTemperature() {
        super();
    }
      
    public final void setAtomColor(Atom a) {
        float red, blue;
        double ke =((AtomType.Wall)a.type).getTemperature();
        if(ke > KEMax) {red = 1.0f;}
        else if(ke < KEMin) {red = 0.0f;}
        else {red = (float)((ke-KEMin)*range);}
        blue = 1.0f - red;
        a.setColor(new Color(red, 0.0f, blue));
    }
}
