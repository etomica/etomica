/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

//class includes a main method to demonstrate and test its use
package etomica.graphics;

import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;
import etomica.data.types.DataDouble;
import etomica.units.dimensions.Null;
import etomica.units.systems.UnitSystem;
import etomica.util.Constants;
import etomica.util.EnumeratedType;

import javax.swing.*;

/**
 * A simple display of a single value in a textbox with an associated label.
 * Value is obtained from an associated DataSource.
 * A label and unit is associated with the value.
 *
 * @author David Kofke
 */
 
public class DisplayTextBox extends Display implements IDataSink, javax.swing.event.ChangeListener {
    
    /**
     * Descriptive text label to be displayed with the value
     */
    protected JLabel jLabel;
    private Constants.CompassDirection labelPosition = Constants.CompassDirection.NORTH;
    /**
     * Object for displaying the value as a text field
     */
    protected JTextField value;
    /**
     * Displayed panel that holds the label and value
     * (not yet used; meant to implement to make lightweight display)
     */
    protected JPanel panel = new JPanel(new java.awt.BorderLayout());
    /**
     * Integer specifying the number of significant figures to be displayed.
     * Default is 4.
     */
    int precision;
    private LabelType labelType;
    private boolean integerDisplay;
    
    /**
     * Physical units associated with the displayed value.
     * Default is null (dimensionless).
     */
    protected etomica.units.Unit unit;
    
    public DisplayTextBox() {
        super();
        this.unit = Null.UNIT;
        jLabel = new JLabel();
        value = new JTextField("");
        value.setEditable(false);
        panel.add(value, java.awt.BorderLayout.CENTER);
        setLabelType(LabelType.STRING);
        setLabel("");
 //       panel.setMinimumSize(new java.awt.Dimension(80,60));
        setPrecision(4);
        setIntegerDisplay(false);
        
 /*       addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent evt) {
                if((evt.getModifiers() & InputEvent.BUTTON3_MASK) != 0) {//right-click
                    Device editor = new DeviceUnitEditor(DisplayBox.this);
                    editor.graphic(null).setLocation(evt.getX(), evt.getY());
                    parentSimulation().add(editor.graphic(null));
                    parentSimulation().validate();
                    parentSimulation().repaint();
                }
            }
        });  */
    }//end of constructor
    
    public void putDataInfo(IDataInfo dataInfo) {
        if (dataInfo.getLength() != 1) {
            throw new RuntimeException("DisplayTextBox wants only 1 value");
        }
        if(unit == Null.UNIT) {
            unit = dataInfo.getDimension().getUnit(UnitSystem.SIM);
        }
        if(label.equals("")) {
            setLabel(dataInfo.getLabel());
        }
    }

    /**
     * calls doUpdate method.  Implementation of ChangeListener interface.
     */
    public void stateChanged(javax.swing.event.ChangeEvent evt) {
        panel.repaint();
    }
    
    /**
     * Accessor method to set the physical units of the displayed value.
     * Text describing unit is used in label.
     */
    public void setUnit(etomica.units.Unit u) {
        unit = u;
    }
    /**
     * Returns the physical units of the displayed value.
     */
    public etomica.units.Unit getUnit() {return unit;}
    
    public java.awt.Component graphic(Object obj) {return panel;}
    
    
    /**
     * @return Returns the integerDisplay.
     */
    public boolean isIntegerDisplay() {
        return integerDisplay;
    }
    /**
     * @param integerDisplay The integerDisplay to set.
     */
    public void setIntegerDisplay(boolean integerDisplay) {
        this.integerDisplay = integerDisplay;
    }
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public int getPrecision() {return precision;}
    
    /**
     * Accessor method of the precision, which specifies the number of significant figures to be displayed.
     */
    public void setPrecision(int n) {
        // the actual size of the box is closer to n+3, which is what we actually want
        value.setColumns(n+2);
        precision = n;
    }
    
    /**
     * Sets the value of a descriptive label using the given string.
     */
    public void setLabel(String s) {
        String suffix = (unit.symbol().length() > 0) ? " ("+unit.symbol()+")" : "";
        super.setLabel(s+suffix);
        jLabel.setText(s+suffix);
        if(labelType == LabelType.BORDER) {
            panel.setBorder(new javax.swing.border.TitledBorder(s+suffix));
        }
        if(labelType == LabelType.STRING) setLabelPosition(labelPosition);
    }
    /**
     * @return the current value of the descriptive label.
     */
    public String getLabel() {return jLabel.getText();}
    
    /**
     * Sets label to the given value if it was not previously set.
     * If setLabel was previously called, this method has no effect.
     * This method is usually invoked automatically when this data
     * sink is attached to a data pipe.
     */
    public void setDefaultLabel(String defaultLabel) {
        if(getLabel() == "") setLabel(defaultLabel);
    }


    public void setLabelType(LabelType labelType) {
        this.labelType = labelType;
        if(labelType != LabelType.BORDER) panel.setBorder(new javax.swing.border.EmptyBorder(2,2,2,2));
        if(labelType != LabelType.STRING) panel.remove(jLabel);
        setLabel(jLabel.getText());
    }
    public LabelType getLabelType() {
        return labelType;
    }

    public void setLabelPosition(Constants.CompassDirection position) {
        labelPosition = position;
        if(labelType != LabelType.STRING) return;
        panel.remove(jLabel);
        panel.add(jLabel,position.toString());//toString() returns the corresponding BorderLayout constant
//        support.firePropertyChange("label",oldLabel,label);
        panel.revalidate();
        panel.repaint();
    }
    
    public Constants.CompassDirection getLabelPosition() {return labelPosition;}
    /**
     * Sets the display text to reflect the desired value from the datasource.
     */
    public void putData(IData data) {
        double xValue = unit.fromSim(((DataDouble)data).x);
        if(integerDisplay) {
            value.setText(Integer.toString((int)xValue));
        } else {
            value.setText(format(xValue,precision));
        }
        panel.repaint();
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {

   //     etomica.simulation.prototypes.HSMD2D sim = new etomica.simulation.prototypes.HSMD2D();
   //     Simulation.instance = sim;

        Simulation sim = Simulation.instance;
	    IntegratorHard integrator = new IntegratorHard();
	    Species species = new SpeciesSpheres();
	    species.setNMolecules(25);
	    Box box = new Box();
	    Controller controller = new Controller();
	    Display display = new DisplayBox();
	    IntegratorMD.Timer timer = integrator.new Timer(integrator.chronoMeter());
	    timer.setUpdateInterval(10);
		
        Potential2 potential = new P2SquareWell(sim);
 //       Potential2 potential = new P2HardSphere(sim);
		sim.elementCoordinator.go();
        Potential2.Agent potentialAgent = (Potential2.Agent)potential.getAgent(box);
        potentialAgent.setIterator(new AtomPairIterator(box));

        //part that is unique to this demonstration
        Modifier mod1 = new Modifier(integrator, "timeStep");
        Meter ke = new MeterKineticEnergy();
        Meter energy = new MeterEnergy();
        ke.setBox(box);
        DisplayBox box = new DisplayBox();
        box.setDatumSource(ke);
        box.setUpdateInterval(10);
        //DisplayBox showing the current value (default is most recent, but this is zero because meter is inactive (not keeping averages), and thus doesn't hold a most-recent value)
        DisplayBox box0 = new DisplayBox(energy);
        box0.setWhichValue(MeterAbstract.CURRENT);
        box0.setLabelType(DisplayBox.STRING);
        box0.setLabelPosition(Constants.SOUTH);
        //here's a DisplayBox tied to a Modifier
		DisplayBox box1 = new DisplayBox();
		box1.setDatumSource(mod1);
        //end of unique part
                                            
		Simulation.instance.elementCoordinator.go(); 		                                    
        Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main  
    */

    
    
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

        if (Double.isNaN(d) || Double.isInfinite(d))
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

        // original = d * 10^exponent
        // 1 < d <= 10
        
        int p = precision;
        while (--p > 0)
            d *= 10;
        
        long ld = Math.round(d);
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
        if (Math.abs(exponent) < precision) {
            while (exponent > 0) {
                ++decimalPoint;
                --exponent;
            }
            while (exponent < 0) {
                --decimalPoint;
                ++exponent;
            }
            // exponent is now 0

            if (decimalPoint < 0) {
                buffer.append("0.");
                while (decimalPoint < -1) {
                    buffer.append("0");
                    ++decimalPoint;
                }
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

    /**
     * Typed constant used to indicate the type of label to be used with the display.
     */
     
	public static class LabelType extends EnumeratedType {
        public LabelType(String label) {super(label);}       
        public static final LabelType BORDER = new LabelType("Border");
        public static final LabelType STRING = new LabelType("String");

        public static final LabelType[] choices() { 
            return new LabelType[] {BORDER,STRING};
        }
    }
    
}
