//class includes a main method to demonstrate and test its use
package etomica;

import etomica.units.*;
import javax.swing.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.InputEvent;

/**
 * A simple display of a single value in a textbox with an associated label.
 * Value is obtained from an associated DatumSource.
 * A label and unit is associated with the value.
 */
 
public class DisplayBox extends Display implements Dimensioned, DatumSource.User, EtomicaElement {
    
    public String getVersion() {return "DisplayBox:01.06.02/"+Display.VERSION;}
    /**
     * Descriptive text label to be displayed with the value
     */
    protected JLabel label;
    private Constants.Direction labelPosition = Constants.NORTH;
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
     * Datum source that generates the displayed value
     */
    protected DatumSource source;
    /**
     * Integer specifying the number of significant figures to be displayed.
     * Default is 4.
     */
    int precision;
    private DataSource.ValueType whichValue;
    private LabelType labelType;
    
    /**
     * Physical units associated with the displayed value.
     * Default is null (dimensionless).
     */
    protected Unit unit;
    
    public DisplayBox() {
        this(Simulation.instance);
    }
    public DisplayBox(DatumSource m) {
        this(Simulation.instance, m);
    }
    public DisplayBox(Simulation sim, DatumSource m) {
        this(sim);
        setDatumSource(m);
    }
    public DisplayBox(Simulation sim) {
        super(sim);
        label = new JLabel("Label");
        value = new JTextField("");
        value.setEditable(false);
        panel.add(value, java.awt.BorderLayout.CENTER);
        setLabelType(STRING);
 //       panel.setMinimumSize(new java.awt.Dimension(80,60));
        unit = new Unit(BaseUnit.Null.UNIT);
        setUpdateInterval(5);
        setPrecision(4);
        
/*        addMouseListener(new MouseAdapter() {
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
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple display of one meter's value with a label");
        return info;
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
        if(source != null) return source.getDimension();
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
    public void setDatumSource(DatumSource m) {
        source = m;
        if(m instanceof Meter && 
            !(whichValue instanceof MeterAbstract.ValueType)) setWhichValue(MeterAbstract.MOST_RECENT);
        setUnit(m.getDimension().defaultIOUnit());
        setLabel();
    }
    
    /**
     * Accessor method for the meter that generates the displayed value.
     */
    public DatumSource getDatumSource() {
        return source;
    }
    
    /**
     * Sets the value of a descriptive label using the meter's label and the unit's symbol (abbreviation).
     */
    private void setLabel() {
        if(source == null) return;
        String suffix = (unit.symbol().length() > 0) ? " ("+unit.symbol()+")" : "";
        setLabel(source.getLabel()+suffix);
    }
    
    /**
     * Sets the value of a descriptive label using the given string.
     */
    public void setLabel(String s) {
        label.setText(s);
        if(labelType == BORDER) {
            panel.setBorder(new javax.swing.border.TitledBorder(s));
        }
        if(labelType == STRING) setLabelPosition(labelPosition);
/*        JLabel oldLabel = label;
        label = new JLabel(s);
        panel.remove(oldLabel);
        panel.add(label, 0);
        support.firePropertyChange("label",oldLabel,label);*/
    }
    /**
     * @return the current value of the descriptive label.
     */
    public String getLabel() {return label.getText();}
    

    public void setLabelType(LabelType labelType) {
        this.labelType = labelType;
        if(labelType != BORDER) panel.setBorder(new javax.swing.border.EmptyBorder(2,2,2,2));
        if(labelType != STRING) panel.remove(label);
        setLabel(label.getText());
    }
    public LabelType getLabelType() {
        return labelType;
    }

    public void setLabelPosition(Constants.Direction position) {
        labelPosition = position;
        if(labelType != STRING) return;
        panel.remove(label);
        panel.add(label,position.toString());//toString() returns the corresponding BorderLayout constant
//        support.firePropertyChange("label",oldLabel,label);
        panel.revalidate();
        panel.repaint();
    }
    
    public Constants.Direction getLabelPosition() {return labelPosition;}
    /**
     * Sets the display text to reflect the desired value from the meter.
     */
    public void doUpdate() {
        if(source == null) return;
        value.setText(format(unit.fromSim(source.value(whichValue)),precision));
    }
    /**
     *   @deprecated  Use setWhichValue instead
     */
    public void setUseCurrentValue(boolean b) {
        if(b) setWhichValue(MeterAbstract.MOST_RECENT);
        else  setWhichValue(MeterAbstract.AVERAGE);
    }
    /**
     * Sets whether meter displays average, current value, last block average, etc.
     */
    public void setWhichValue(DataSource.ValueType type) {
        whichValue = type;
        if(source == null) return;
//        if(!(whichValue==MeterAbstract.ValueType.MOST_RECENT) && !meter.isActive())
//            {System.out.println("Warning: setting to use averages but meter is not active");}
    }
    public DataSource.ValueType getWhichValue() {return whichValue;}
    
    public void setMeter(Meter m) {
        setDatumSource(m);
    }
    public Meter getMeter() {return (Meter)getDatumSource();}
    
    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);

        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;

        //part that is unique to this demonstration
        PotentialSquareWell potential = new PotentialSquareWell(sim);
        sim.p2 = new P2SimpleWrapper(sim,sim.potential);
        Modulator mod1 = new Modulator(sim.integrator, "timeStep");
        Meter ke = new MeterKineticEnergy();
        ke.setPhase(sim.phase);
        DisplayBox box = new DisplayBox();
        box.setDatumSource(ke);
        box.setUpdateInterval(10);
        //DisplayBox showing the current value (default is most recent, but this is zero because meter is inactive (not keeping averages), and thus doesn't hold a most-recent value)
        DisplayBox box0 = new DisplayBox.Energy(sim.phase);
        box0.setWhichValue(MeterAbstract.CURRENT);
        box0.setLabelType(DisplayBox.STRING);
        box0.setLabelPosition(Constants.SOUTH);
        //here's a DisplayBox tied to a Modulator
		DisplayBox box1 = new DisplayBox();
		box1.setDatumSource(mod1);
        //end of unique part
                                            
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        f.getContentPane().add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }

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

    /**
     * Typed constant used to indicate the type of label to be used with the display.
     */
     
	public static class LabelType extends Constants.TypedConstant {
        public LabelType(String label) {super(label);}       
        public Constants.TypedConstant[] choices() {return (Constants.TypedConstant[])CHOICES;}
        public static final LabelType[] CHOICES = 
            new LabelType[] {
                new LabelType("Border"),
                new LabelType("String")};
    }//end of LabelType
    public static final LabelType BORDER = LabelType.CHOICES[0];
    public static final LabelType STRING = LabelType.CHOICES[1];
    
}