//class includes a main method to demonstrate and test its use
package etomica;

import etomica.units.*;
import javax.swing.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.InputEvent;

/**
 * A simple display of a single value in a textbox with an associated label.
 * Value is obtained from an associated Meter, and may be an instantaneous or average value.
 * A label and unit is associated with the value.
 */
 
public class DisplayBox extends Display implements Dimensioned, Meter.User {
    
    /**
     * Descriptive text label to be displayed with the value
     */
    protected JLabel label;
    /**
     * Object for displaying the value as a text field
     */
    protected JTextField value;
    /**
     * Displayed panel that holds the label and value
     * (not yet used; meant to implement to make lightweight display)
     */
    protected JPanel panel = new JPanel();
    /**
     * Meter that generates the displayed value
     */
    protected Meter meter;
    /**
     * Integer specifying the number of significant figures to be displayed.
     * Default is 4.
     */
    int precision;
    private boolean useCurrentValue;  //need to set up flags for recent, average, current
    /**
     * Physical units associated with the displayed value.
     * Default is null (dimensionless).
     */
    protected Unit unit;
    
    public DisplayBox() {
        this(Simulation.instance);
    }
    public DisplayBox(Simulation sim) {
        super(sim);
        label = new JLabel("Label");
        value = new JTextField("");
        value.setEditable(false);
        setPrecision(4);
        panel.add(label);
        panel.add(value);
        unit = new Unit(BaseUnit.Null.UNIT);
        setUseCurrentValue(true);
        setUpdateInterval(5);
        panel.setSize(100,60);
        
        addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent evt) {
                if((evt.getModifiers() & InputEvent.BUTTON3_MASK) != 0) {//right-click
                    Device editor = new DeviceUnitEditor(DisplayBox.this);
                    editor.graphic(null).setLocation(evt.getX(), evt.getY());
                    parentSimulation().add(editor.graphic(null));
                    parentSimulation().validate();
                    parentSimulation().repaint();
                }
            }
        });
        
    }
    
    /**
     * Accessor method to set the physical units of the displayed value.
     * Text describing unit is used in label.
     */
    public void setUnit(Unit u) {
        unit = u;
        setLabel();
    }
    /**
     * Returns the physical units of the displayed value.
     */
    public Unit getUnit() {return unit;}
    
    /**
     * Returns the dimensions of the quantity being measured.
     * Obtained from the meter associated with this display.
     */
    public Dimension dimension() {
        if(meter != null) return meter.getDimension();
        else return Dimension.NULL;
    }
    
    public java.awt.Component graphic(Object obj) {return panel;}
    
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public int getPrecision() {return precision;}
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public void setPrecision(int n) {
        value.setColumns(n);
        precision = n;
    }
    
    /**
     * Specifies the meter that generates the displayed value.
     */
    public void setMeter(Meter m) {
        meter = m;
        setUnit(m.defaultIOUnit());
        setLabel();
    }
    
    /**
     * Accessor method for the meter that generates the displayed value.
     */
    public Meter getMeter() {
        return meter;
    }
    
    /**
     * Sets the value of a descriptive label using the meter's label and the unit's symbol (abbreviation).
     */
    private void setLabel() {
        if(meter == null) return;
        String suffix = (unit.symbol().length() > 0) ? " ("+unit.symbol()+")" : "";
        setLabel(meter.getLabel()+suffix);
    }
    
    /**
     * Sets the value of a descriptive label using the given string.
     */
    public void setLabel(String s) {
        JLabel oldLabel = label;
        label = new JLabel(s);
        panel.remove(oldLabel);
        panel.add(label, 0);
        support.firePropertyChange("label",oldLabel,label);
    }
    /**
     * @return the current value of the descriptive label.
     */
    public String getLabel() {return label.getText();}
    
    /**
     * Sets the display text to reflect the current or average value from the meter.
     */
    public void doUpdate() {
        if(meter == null) return;
        if(useCurrentValue) {
            value.setText(format(unit.fromSim(meter.mostRecent()),precision));
        }
        else {
            value.setText(format(unit.fromSim(meter.average()),precision));
        }
    }
    
    /**
     * Sets flag indicating if plot should be of instantaneous (current) value or running average
     */
    public void setUseCurrentValue(boolean b) {
        useCurrentValue = b;
        if(meter == null) return;
        if(!useCurrentValue && !meter.isActive()) {System.out.println("Warning: setting to use averages but meter is not active");}
    }
    /**
     * Accessor method for currentValue field, which specifies if display should present current value
     * or running average.
     */
    public boolean getUseCurrentValue() {return useCurrentValue;}

    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        Simulation.makeSimpleSimulation();  //for more general simulations, replace this call with
                                            //construction of the desired pieces of the simulation
        //part that is unique to this demonstration
        Meter ke = new MeterKineticEnergy();
        Phase phase = Simulation.instance.phase(0);
        ke.setPhase(phase);
        DisplayBox box = new DisplayBox();
        box.setMeter(ke);
        box.setUpdateInterval(10);
        //another display box; note that we don't need to save instance to have it show in simulation,
        //but we can't change its properties since we don't keep a handle here
        new DisplayBox.Energy(phase);  
        //end of unique part
                                            
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }

    
    //*******  Some convenience subclasses of DisplayBox  ************//
    
    /**
     * A display box that presents the current value of the total energy of a phase
     */
    public static class Energy extends DisplayBox {
    
        public Energy(Phase phase) {
            super();
            this.setMeter(phase.energy);
//            this.setUnit(Simulation.unitSystem().energy());
        }
    }
    
    
    /********** Utility method for formatting a double to a string **************/
    
    /**
     * Formats a double with a specified number of digits.
     * When java converts a <tt>double</tt> to a <tt>String</tt>
     * it retains the full precision of the number. This can
     * generate 15 decimal places! This method truncates this output
     * to some specified number of decimal places.
     * @param d the double to format
     * @param precision the number of digits desired
     * @return returns the formatted string
     *
     * Taken from the comphys package of Richard Gonsalves of the
     * SUNY Buffalo Department of Physics
     */

    public static String format (double d, int precision) {

        if (d == Double.NaN ||
            d == Double.POSITIVE_INFINITY ||
            d == Double.NEGATIVE_INFINITY)
            return Double.toString(d);
        
        StringBuffer buffer = new StringBuffer(20);
        
        if (d < 0) {
            d = -d;
            buffer.append('-');
        }

        if (d == 0) {
            buffer.append("0.0");
            for (int p = 0; p < precision - 1; p++)
                buffer.append('0');
            return buffer.toString();
        }

        int exponent = 0;
        while (d >= 10) {
            ++exponent;
            d /= 10;
        }
        while (d < 1) {
            --exponent;
            d *= 10;
        }

        if (precision < 0)
            precision = -precision;
        int p = precision;
        while (--p > 0)
            d *= 10;
        long ld = (long) Math.round(d);
        char[] digits = new char[precision];
        p = precision;
	long ld_div_10 = 0;
	long ld_save = ld;
        while (--p >= 0) {
	    ld_div_10 = ld / 10;
            digits[p] = (char) ('0' + ( ld - (ld_div_10 * 10) ));
            ld = ld_div_10;
        }
	if (ld_div_10 > 0) {
	    ld = ld_save / 10;
	    p = precision;
	    while (--p >= 0) {
		ld_div_10 = ld / 10;
		digits[p] = (char) ('0' + ( ld - (ld_div_10 * 10) ));
		ld = ld_div_10;
	    }
	    ++exponent;
	}

        int decimalPoint = 0;
        if (Math.abs(exponent) < 6 || Math.abs(exponent) < precision) {
            while (exponent > 0) {
                ++decimalPoint;
                --exponent;
            }
            while (exponent < 0) {
                --decimalPoint;
                ++exponent;
            }
        }

        if (decimalPoint < 0) {
            buffer.append("0.");
            while (decimalPoint < -1) {
                buffer.append("0");
                ++decimalPoint;
            }
        }

        for (p = 0; p < precision; p++) {
            buffer.append(digits[p]);
            if (p == decimalPoint)
                if (p < precision - 1)
                    buffer.append(".");
        }

        if (exponent != 0)
            buffer.append("E" + exponent);

        return buffer.toString();

    }

}