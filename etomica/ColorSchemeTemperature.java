package simulate;
import java.awt.Color;
/**
 *
 * @author David Kofke
 *
 */
 
public class ColorSchemeTemperature extends ColorScheme {
    
    double TLow, THigh;
    protected double KEMin, KEMax, range;
    
    public ColorSchemeTemperature() {
        this(200.0, 400.0);
    }
    public ColorSchemeTemperature(double TLow, double THigh) {
        setTLow(TLow);
        setTHigh(THigh);
    }
      
    public double getTLow() {return TLow;}
    public void setTLow(double t) {
        TLow = t;
        KEMin = t/Constants.KE2T;
        range = 1.0/(KEMax-KEMin);
    }
    public double getTHigh() {return THigh;}
    public void setTHigh(double t) {
        THigh = t;
        KEMax = t/Constants.KE2T;
        range = 1.0/(KEMax-KEMin);
    }
        
    public final void setAtomColor(Atom a) {
        float red, blue;
        double ke = a.kineticEnergy();
        if(ke > KEMax) {blue = 0.0f;}
        else if(ke < KEMin) {blue = 1.0f;}
        else {blue = (float)((ke-KEMin)*range);}
        red = 1.0f - blue;
        a.setColor(new Color(red, 0.0f, blue));
    }
}
