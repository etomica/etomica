//This class includes a main method to demonstrate its use
package simulate;
import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.InputEvent;
import simulate.units.*;
import simulate.utility.StringUtility;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

/**
 * Device the changes a property using a graphical slider, via a Modulator.
 * Entire visible device consists of slider and a label.  The slider is an instance
 * of a JSlider class, and the label is a JLabel instance; the properties of either object
 * may be set directly by accessing the objects with the getSlider() or getLabel() methods.
 *
 * @see Modulator
 */
public class DeviceSlider extends Device {
    
    /**
     * Descriptive text label to be displayed with the value
     */
    protected JLabel label;
    /**
     * Modulator connecting the slider to the property
     */
    protected Modulator modulator;
    /**
     * Swing slider displayed to screen
     */
    protected JSlider slider;
    /**
     * Panel to hold slider and label, and returned by the graphic method
     */
    protected JPanel panel;
    /**
     * Object with property being modulated
     */
    protected Object component;
    /**
     * Property being modulated
     */
    protected String property;
    
    private int minimum, maximum;
    
    public DeviceSlider() {
        this(Simulation.instance);
    }
    public DeviceSlider(Simulation sim) {
        super(sim);
        slider = new JSlider();
        slider.setSize(200,40);
        slider.setPaintTicks(true);
        slider.setPaintLabels(true);
        slider.setValue(300);
        slider.setMinimum(100);
        slider.setMaximum(500);
        slider.setMajorTickSpacing(100);
//        slider.setMinorTickSpacing(50);
        slider.addChangeListener(new SliderListener());  //SliderListener is an inner class defined below
        label = new JLabel("");
        label.setHorizontalAlignment(SwingConstants.CENTER);
        panel = new JPanel();
        panel.setLayout(new BorderLayout(0,0));
        panel.add(slider, BorderLayout.NORTH);        
        panel.add(label, BorderLayout.CENTER);
        
        label.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent evt) {
                if((evt.getModifiers() & InputEvent.BUTTON3_MASK) != 0) {//right-click
                    Device editor = new DeviceUnitEditor(DeviceSlider.this);
                 //   editor.graphic(null).setLocation(evt.getX(), evt.getY());
                 //   panel.getParent().add(editor.graphic(null));
                    parentSimulation().add(editor.graphic(null));
                    parentSimulation().validate();
                    parentSimulation().repaint();
                }
            }
        });
        
    }
    
    /**
     * Override superclass setUnit method to update label when unit is changed
     */
    public void setUnit(Unit u) {
        super.setUnit(u);
        setLabelDefault();
    }
    
    /**
     * Constructs a slider connected to the given property of the given object
     */
    public DeviceSlider(Object object, String property) {
        this(new Modulator(object, property));
        component = object;
        this.property = property;
    }
    /**
     * Constructs a slider connected to the get/set Value methods of the given Modulator
     */
    public DeviceSlider(Modulator m) {
        this();
        //set component and property in some way
        setModulator(m);
    }
    
    public final void setModulator(Modulator m) {
        modulator = m;
        unit = modulator.getDimension().defaultIOUnit();
        setLabelDefault();
        slider.setValue((int)unit.fromSim(modulator.getValue()));
        setMinimum(getMinimum());
        setMaximum(getMaximum());
    }
    public final Modulator getModulator() {return modulator;}
    
    public String getProperty() {return property;}
    public void setProperty(String s) {
        property = s;
        if(component != null) setModulator(new Modulator(component, property));
    }
    
    public Object getComponent() {return component;}
    public void setComponent(Object obj) {
        component = obj;
        if(property != null) setModulator(new Modulator(component,property));
    }
    
    public void setMinimum(int min) {
        minimum = min;
        slider.setMinimum(min);
        slider.setInverted(maximum < minimum);
        setTicks();
    }
    public int getMinimum() {return minimum;}
    
    public void setMaximum(int max) {
        maximum = max;
        slider.setMaximum(max);
        slider.setInverted(maximum < minimum);
        setTicks();
    }
    public int getMaximum() {return maximum;}
    
    private void setTicks() {
        int spacing = (getMaximum()-getMinimum())/4;
        if(spacing <= 0) return;
        slider.setMajorTickSpacing(spacing);
        slider.setMinorTickSpacing(slider.getMajorTickSpacing()/2);
        //need to do the following because JSlider does not automatically
        //reset labels if they have been set before
        slider.setLabelTable(slider.createStandardLabels(spacing));
    }
    
    /**
     * Returns the GUI element for display in the simulation.
     * Consists of a panel containing the swing slider and label.
     */
    public java.awt.Component graphic(Object obj) {
        return panel;
    }
    
    
    /**
     * Sets the value of a descriptive label using the meter's label and the unit's symbol (abbreviation).
     */
    private void setLabelDefault() {
        String suffix = (unit.symbol().length() > 0) ? " ("+unit.symbol()+")" : "";
        if(modulator != null) 
            label.setText(StringUtility.capitalize(modulator.getProperty())+suffix);
    }

    /**
     * Sets the value of a descriptive label using the given string.
     */
    public void setLabel(JLabel l) {label = l;}
    /**
     * @return the current instance of the descriptive label.
     */
    public JLabel getLabel() {return label;}
    
    /**
     * @return a handle to the JSlider instance used by this slider device
     */
    public JSlider getSlider() {return slider;}
    
    /**
     * Performs no action (slider can not be accessed for setting).
     * Method exists so that slider and its properties (which can be set) will
     * show up on the property sheet.
     */
     public void setSlider(JSlider s) {}
    
    /**
     * Slider listener, which relays the slider change events to the modulator
     */
    private class SliderListener implements ChangeListener {
        public void stateChanged(ChangeEvent evt) {
            if(modulator!=null) modulator.setValue(unit.toSim(slider.getValue()));
       }
    }
    
    /**
     * Method to demonstrate and test the use of this class.  
     * Slider is used to control the temperature of a hard-disk MD simulation
     */
    public static void main(String[] args) {
        Frame f = new Frame();   //create a window
        f.setSize(600,350);
        Simulation.makeSimpleSimulation();  
        
        //here's the part unique to this class
        Integrator integrator = Simulation.integrator(0);
        DeviceSlider mySlider = new DeviceSlider(integrator,"temperature");
        integrator.setIsothermal(true);
        //end of unique part
 
        Simulation.instance.elementCoordinator.go();
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(WindowEvent e) {System.exit(0);}
        });
    }
}