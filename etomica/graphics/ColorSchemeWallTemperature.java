package etomica.graphics;
import java.awt.Color;

import etomica.Atom;
import etomica.atom.AtomTypeWall;
/**
 * Sets color of a wall based on the value of its temperature field.
 * @author David Kofke
 *
 */
 
public class ColorSchemeWallTemperature extends ColorSchemeTemperature {
    
    public ColorSchemeWallTemperature() {
        super();
    }
      
    public final Color atomColor(Atom a) {
        float red, blue;
        double ke =((AtomTypeWall)a.type).getTemperature();
        if(ke > KEMax) {red = 1.0f;}
        else if(ke < KEMin) {red = 0.0f;}
        else {red = (float)((ke-KEMin)*range);}
        blue = 1.0f - red;
        return new Color(red, 0.0f, blue);
    }
}
