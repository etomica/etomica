//class includes a main method to demonstrate and test its use
package etomica.graphics;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import etomica.Constants;
import etomica.DataSource;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;

/**
 * A simple display of a single value in a textbox with an associated label.
 * Value is obtained from an associated DataSource.
 * A label and unit is associated with the value.
 *
 * @author David Kofke
 */
 
public class DisplayBox extends Display implements etomica.units.Dimensioned, DataSource.User, EtomicaElement, javax.swing.event.ChangeListener {
    
    /**
     * Descriptive text label to be displayed with the value
     */
    protected JLabel jLabel;
    private Constants.CompassDirection labelPosition = Constants.NORTH;
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
    protected DataSource source;
    /**
     * Integer specifying the number of significant figures to be displayed.
     * Default is 4.
     */
    int precision;
    private LabelType labelType;
    
    /**
     * Physical units associated with the displayed value.
     * Default is null (dimensionless).
     */
    protected etomica.units.Unit unit;
    
    public DisplayBox(DataSource m) {
        this();
        setDataSource(m);
    }
    public DisplayBox() {
        super();
        jLabel = new JLabel("Label");
        value = new JTextField("");
        value.setEditable(false);
        panel.add(value, java.awt.BorderLayout.CENTER);
        setLabelType(STRING);
 //       panel.setMinimumSize(new java.awt.Dimension(80,60));
        unit = new etomica.units.PrefixedUnit(etomica.units.BaseUnit.Null.UNIT);
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
        EtomicaInfo info = new EtomicaInfo("Simple display of one data source's value with a label");
        return info;
    }
    
    /** 
     * calls doUpdate method.  Implementation of ChangeListener interface.
     */
    public void stateChanged(javax.swing.event.ChangeEvent evt) {
        doUpdate();
    }
    
    /**
     * Accessor method to set the physical units of the displayed value.
     * Text describing unit is used in label.
     */
    public void setUnit(etomica.units.Unit u) {
        unit = u;
        setLabel();
    }
    /**
     * Returns the physical units of the displayed value.
     */
    public etomica.units.Unit getUnit() {return unit;}
    
    /**
     * Returns the dimensions of the quantity being measured.
     * Obtained from the datsource associated with this display.
     */
    public etomica.units.Dimension dimension() {
        if(source != null) return source.getDimension();
        return etomica.units.Dimension.NULL;
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
     * Specifies the datasourcethat generates the displayed value.
     */
    public void setDataSource(DataSource m) {
        source = m;
        if (unit==null) {
            setUnit(m.getDimension().defaultIOUnit());
        }
        setLabel();
    }
    
    /**
     * Accessor method for the datsource that generates the displayed value.
     */
    public DataSource getDataSource() {
        return source;
    }
    
    /**
     * Sets the value of a descriptive label using the datsource's label and the unit's symbol (abbreviation).
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
        jLabel.setText(s);
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
    public String getLabel() {return jLabel.getText();}
    

    public void setLabelType(LabelType labelType) {
        this.labelType = labelType;
        if(labelType != BORDER) panel.setBorder(new javax.swing.border.EmptyBorder(2,2,2,2));
        if(labelType != STRING) panel.remove(jLabel);
        setLabel(jLabel.getText());
    }
    public LabelType getLabelType() {
        return labelType;
    }

    public void setLabelPosition(Constants.CompassDirection position) {
        labelPosition = position;
        if(labelType != STRING) return;
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
    public void doUpdate() {
        if(source == null) return;
        value.setText(format(unit.fromSim(source.getData()[0]),precision));
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {

   //     etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
   //     Simulation.instance = sim;

        Simulation sim = Simulation.instance;
	    IntegratorHard integrator = new IntegratorHard();
	    Species species = new SpeciesSpheres();
	    species.setNMolecules(25);
	    Phase phase = new Phase();
	    Controller controller = new Controller();
	    Display display = new DisplayPhase();
	    IntegratorMD.Timer timer = integrator.new Timer(integrator.chronoMeter());
	    timer.setUpdateInterval(10);
		
        Potential2 potential = new P2SquareWell(sim);
 //       Potential2 potential = new P2HardSphere(sim);
		sim.elementCoordinator.go();
        Potential2.Agent potentialAgent = (Potential2.Agent)potential.getAgent(phase);
        potentialAgent.setIterator(new AtomPairIterator(phase));

        //part that is unique to this demonstration
        Modifier mod1 = new Modifier(integrator, "timeStep");
        Meter ke = new MeterKineticEnergy();
        Meter energy = new MeterEnergy();
        ke.setPhase(phase);
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
        public Constants.TypedConstant[] choices() {return CHOICES;}
        public static final LabelType[] CHOICES = 
            new LabelType[] {
                new LabelType("Border"),
                new LabelType("String")};
    }//end of LabelType
    public static final LabelType BORDER = LabelType.CHOICES[0];
    public static final LabelType STRING = LabelType.CHOICES[1];
    
}