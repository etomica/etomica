//includes a main method to demonstrate use and to test
package etomica.graphics;
import etomica.*;
import etomica.units.*;
import java.awt.Color;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.Frame;
/**
 * Colors atoms according to their kinetic energy.
 * Atoms with high KE are colored red, and those with low KE are colored blue.
 * Range of low..high is adjustable.
 *
 * @author David Kofke
 *
 */
 
public class ColorSchemeTemperature extends ColorScheme {
    
    double TLow, THigh;
    protected double KEMin, KEMax, range;
    
    /**
     * Constructs with default low of 200K and high of 400K.
     */
    public ColorSchemeTemperature() {
        this(Kelvin.UNIT.toSim(200.), Kelvin.UNIT.toSim(400.));
    }
    public ColorSchemeTemperature(double TLow, double THigh) {
        setTLow(TLow);
        setTHigh(THigh);
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
        
    public void colorAtom(Atom a) {
        float red, blue;
        double ke = a.coord.kineticEnergy();
        if(ke > KEMax) {blue = 0.0f;}
        else if(ke < KEMin) {blue = 1.0f;}
        else {blue = (float)((ke-KEMin)*range);}
        red = 1.0f - blue;
        a.setColor(new Color(red, 0.0f, blue));
    }
    
    /**
     * Demonstrates how this class is implemented.  
     * Appropriate for system having only one species.  For multiple species, see ColorSchemeBySpecies
     */
/*    public static void main(String[] args) {
        Frame f = new Frame();   //create a window
        f.setSize(600,350);
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        
        //part unique to this example
             //get handles to components we need
        DisplayPhase display = sim.display;
             //instantiate color scheme and link appropriately
        ColorSchemeTemperature colorScheme = new ColorSchemeTemperature();
        colorScheme.setTLow(Kelvin.UNIT.toSim(100.)); //select low of temperature scale
        
        display.setColorScheme(colorScheme);
        //end of unique part
		Simulation.instance.elementCoordinator.go(); 
        f.add(Simulation.instance.panel());         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(WindowEvent e) {System.exit(0);}
        });
    }
  */  
}
