//includes a main method to demonstrate use and to test
package etomica.graphics2;

import etomica.atom.Atom;
import etomica.atom.AtomTypeLeaf;
import etomica.space.ICoordinateKinetic;
import etomica.units.Dimension;
import etomica.units.Kelvin;
/**
 * Colors atoms according to their kinetic energy.
 * Atoms with high KE are colored red, and those with low KE are colored blue.
 * Range of low..high is adjustable.
 *
 * @author David Kofke
 *
 */
 
public class ColorSchemeTemperature implements ColorScheme {
    
    double TLow, THigh;
    protected double KEMin, KEMax, range;
    protected int num_colors;
    
    /**
     * Constructs with default low of 200K and high of 400K.
     */
    public ColorSchemeTemperature( int ncolors ) {
        this( ncolors, Kelvin.UNIT.toSim(200.), Kelvin.UNIT.toSim(400.));
     
    }
    public ColorSchemeTemperature(  int ncolors, double TLow, double THigh ) {
        setTLow(TLow);
        setTHigh(THigh);
        num_colors = ncolors;
    }
      
    public double getTLow() {return TLow;}
    public void setTLow(double t) {
        TLow = t;
        KEMin = t;
        range = 1.0/(KEMax-KEMin);
    }
    public Dimension getTLowDimension() {return Dimension.TEMPERATURE;}
    public double getTHigh() {return THigh;}
    public void setTHigh(double t) {
        THigh = t;
        KEMax = t;
        range = 1.0/(KEMax-KEMin);
    }
        
    public int atomColor(Atom a) {
        float blueness = 0.0f;
        double ke = ((AtomTypeLeaf)a.type).getMass()*((ICoordinateKinetic)a.coord).velocity().squared();
        if(ke > KEMax) {blueness = 0.0f;}
        else if(ke < KEMin) {blueness = 1.0f;}
        else {blueness = (float)((KEMax-ke)*range);}

        int color = (int)( num_colors*blueness );
        return color;
    }

	public int getNumColors() {
		return num_colors;
	}
	public Color getColor(int index) {
		float blue = ((float)index)/num_colors;
		float red  = 1.0f - blue;
		float green = 0;
		return new Color( red, green, blue);
	}
}
