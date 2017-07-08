/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

//includes a main method to demonstrate use and to test
package etomica.graphics;
import java.awt.Color;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.units.dimensions.Dimension;
import etomica.units.Kelvin;
import etomica.units.dimensions.Temperature;
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
    	super();
        setTLow(TLow);
        setTHigh(THigh);
    }
      
    public double getTLow() {return TLow;}
    public void setTLow(double t) {
        TLow = t;
        KEMin = t;
        range = 1.0/(KEMax-KEMin);
    }
    public Dimension getTLowDimension() {return Temperature.DIMENSION;}
    public double getTHigh() {return THigh;}
    public void setTHigh(double t) {
        THigh = t;
        KEMax = t;
        range = 1.0/(KEMax-KEMin);
    }
        
    public Color getAtomColor(IAtom a) {
        float red, blue;
        double ke = a.getType().getMass()*((IAtomKinetic)a).getVelocity().squared();
        if(ke > KEMax) {blue = 0.0f;}
        else if(ke < KEMin) {blue = 1.0f;}
        else {blue = (float)((KEMax-ke)*range);}
        red = 1.0f - blue;
        return new Color(red, 0.0f, blue);
    }
    
    /**
     * Demonstrates how this class is implemented.  
     * Appropriate for system having only one species.  For multiple species, see ColorSchemeBySpecies
     */
/*    public static void main(String[] args) {
        Frame f = new Frame();   //create a window
        f.setSize(600,350);
        etomica.simulation.prototypes.HSMD2D sim = new etomica.simulation.prototypes.HSMD2D();
        Simulation.instance = sim;
        
        //part unique to this example
             //get handles to components we need
        DisplayBox display = sim.display;
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
