package simulate;
import java.awt.Color;
/**
 *
 * @author David Kofke
 *
 */
 
public class ColorSchemeWallTemperature extends ColorSchemeTemperature {
    
    public ColorSchemeWallTemperature() {
        super();
    }
      
    public final void setAtomColor(AtomWall a) {
        float red, blue;
        double ke =a.getTemperature();
        if(ke > KEMax) {red = 1.0f;}
        else if(ke < KEMin) {red = 0.0f;}
        else {red = (float)((ke-KEMin)*range);}
        blue = 1.0f - red;
        a.setColor(new Color(red, 0.0f, blue));
    }
}
