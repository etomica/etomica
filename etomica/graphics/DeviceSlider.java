//This class includes a main method to demonstrate its use
package etomica.graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.action.Action;
import etomica.action.activity.Controller;
import etomica.event.ChangeEventManager;
import etomica.modifier.Modifier;
import etomica.modifier.ModifierGeneral;
import etomica.modifier.ModifyAction;
import etomica.units.Unit;
import etomica.util.StringUtility;

/**
 * Device the changes a property using a graphical slider, via a Modifier.
 *
 * @see ModifierGeneral
 */
 
 /* History
  * 10/08/02 (SKK) modified for DecimalSlider, textBox
  * 10/12/02 (DAK) added init method
  * 10/13/02 (DAK) restored nMajor and its accessor/mutator methods
  *                changed graphic method to always return panel, never just slider
  *                always adding slider to panel
  */
  
public class DeviceSlider extends Device implements EtomicaElement {
    
    /**
     * Descriptive text label to be displayed with the value
     * No longer used -- instead apply label via a title border on the slider itself
     */
    private String label;
    /**
     * Modifier connecting the slider to the property
     */
    protected ModifyAction modifyAction;
    /**
     * Subclass of Swing slider displayed to screen
     * located in utility package 
     */
    protected DecimalSlider slider;
    /**
     * Object with property being modulated
     */
    protected Object component;
    /**
     * Property being modulated
     */
    protected String property;
    /**
     * Values for maximum and minimum values of slider in double form
     */
    private double minimum, maximum;
    /**
     * boolean showValues for showing values on textfield
     * boolean editvalus for editing the values in textfield, which change slider tip
     */
    private boolean showValues, editValues;
    /*
     * To show slider with value in one 
     */
    private JPanel panel; 
    /**
     * To show the values of slider
     */
    private JTextField textField;
    
    /** 
     * column of textfield to show the value of slider correctly with horizontal view
     * default value is five and affected if precision is greater than 5
     */    
    private int column;
    
    /**
     * Layout instance to show slider and textfield 
     */    
    private GridBagLayout gbLayout;    
    private GridBagConstraints gbConst; 
    private boolean showBorder = false;
    private int nMajor = 3;
    protected Action targetAction;
    
    protected final ChangeEventManager changeEventManager = new ChangeEventManager(this);
    
    public DeviceSlider(Controller controller) {
        super(controller);
        init();
    }
    
    /**
     * Constructs a slider connected to the given property of the given object
     */
    public DeviceSlider(Controller controller, Object object, String property) {
        this(controller, new ModifierGeneral(object, property));
        component = object;
        this.property = property;
    }
    /**
     * Constructs a slider connected to the get/set Value methods of the given Modifier
     */
    public DeviceSlider(Controller controller, Modifier m) {
        this(controller);
        //set component and property in some way
        setModifier(m);
    }

    private void init() {
        textField = new JTextField("");
        textField.setFont(new java.awt.Font("",0,15));
        textField.setHorizontalAlignment(JTextField.CENTER);
        gbLayout = new GridBagLayout();
        gbConst = new GridBagConstraints();    
        panel = new JPanel();  
//        panel.setBorder(new javax.swing.border.TitledBorder("JPST")); //JPanel of Slider and TextField
        setLabel("");
        panel.setLayout(gbLayout);
        column = 5;
        slider = new DecimalSlider();
        slider.setSize(200,40);
        slider.setPaintTicks(true);
        slider.setPaintLabels(true);
        slider.setDecimalSliderValue(300);
        setMinimum(100);
        setMaximum(500);
        slider.setDecimalSliderMajorTickSpacing(100);
        slider.setDecimalSliderMinorTickSpacing(50);

        gbConst.gridx = 0; gbConst.gridy = 0;
        gbLayout.setConstraints(getSlider(), gbConst);
        panel.add(getSlider());        

        setShowValues(false); // default is false to show values of slider
        setEditValues(false); // default is false to edit values of slider thru textField
    }
    
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo();
        info.setDescription("Slider-type device for changing a property");
        return info;
    }

    /**
     * Override superclass setUnit method to update label when unit is changed
     */
    public void setUnit(Unit u) {
        super.setUnit(u);
        setLabelDefault();
    }
    
    public final void setModifier(Modifier m) {
        if(m == null) throw new NullPointerException();
        modifyAction = null;
        if(unit == null) {
            setUnit(m.getDimension().defaultIOUnit());
        }
        slider.setDecimalSliderValue(unit.fromSim(m.getValue()));
        modifyAction = new ModifyAction(m);
        targetAction = modifyAction;//need to keep this distinct from modifyAction, in case subclasses want to do more than modifyAction when slider is moved
        setLabelDefault();
        setMinimum(getMinimum());
        setMaximum(getMaximum());
    }
    public final Modifier getModifier() {return modifyAction.getWrappedModifier();}
    
    public String getProperty() {return property;}
    public void setProperty(String s) {
        property = s;
        if(component != null) setModifier(new ModifierGeneral(component, property));
    }
    
    public Object getComponent() {return component;}
    public void setComponent(Object obj) {
        component = obj;
        if(property != null) setModifier(new ModifierGeneral(component,property));
    }
    
    public double getMinimum() {return minimum;}
    /**
     * Sets minimum value of slider; should be called after
     * any calls to setPrecision.
     */
    public void setMinimum(double min) {
        minimum = min;
        ModifyAction tmpModifier = modifyAction;
        modifyAction = null;
        slider.setDecimalSliderMinimum(min);
        slider.setInverted(maximum < minimum);
        setTicks();
        modifyAction = tmpModifier;
    }
        
    public double getMaximum() {return maximum;}
    /**
     * Sets maximum value of slider; should be called after
     * any calls to setPrecision.
     */
    public void setMaximum(double max) {
        maximum = max;
        ModifyAction tmpModifier = modifyAction;
        modifyAction = null;
        slider.setDecimalSliderMaximum(max);
        slider.setInverted(maximum < minimum);
        setTicks();
        modifyAction = tmpModifier;
    }
    
    private boolean showMinorValues = false;
    public void setShowMinorValues(boolean b){
        showMinorValues = b;
    }
    
    public void setNMajor(int n) {
        nMajor = n;
        setTicks();
    }
    public int getNMajor() {return nMajor;}
    
    private void setTicks() {
        double minorTick = 1.0 ;
       if(showMinorValues){ minorTick = 2;}
       double spacing = (getMaximum()-getMinimum())/nMajor; 
       if(spacing <= 0) return;
        slider.setDecimalSliderMajorTickSpacing(spacing);
        slider.setDecimalSliderMinorTickSpacing(spacing/2.0);
        //need to do the following because JSlider does not automatically
        //reset labels if they have been set before
        slider.setDecimalSliderLabelTable(slider.createDecimalSliderStandardLabels(Math.max(spacing/minorTick,1)));
    }
    
    public boolean getShowValues(){ return showValues;}
    public void setShowValues(boolean b){
        showValues = b;
        if(showValues){ 
             setSliderValueShape("VERTICAL");
             textField.addActionListener(new ActionListener(){
                  public void actionPerformed( ActionEvent e){
                    //sets value of slider, which then fires event to modifier, taking care of units
                       setValue((Double.parseDouble(e.getActionCommand())));
//                       setValue(unit.toSim(Double.parseDouble(e.getActionCommand())));
             }});}     
    }

    public void setSliderValueShape(String s){
        panel.removeAll();
        textField.setColumns(column);                                                                  
        gbConst.gridx = 0; gbConst.gridy = 0;
        gbLayout.setConstraints(getSlider(), gbConst);
        panel.add(getSlider());        
        if(s=="HORIZONTAL") {
            gbConst.gridx = 1; gbConst.gridy = 0; 
            gbLayout.setConstraints(textField, gbConst);
            panel.add(textField);}
        if(s=="VERTICAL") { 
            gbConst.fill = GridBagConstraints.HORIZONTAL;
            gbConst.gridx = 0; gbConst.gridy = 1; 
            gbLayout.setConstraints(textField, gbConst);
            panel.add(textField);
        }
    }
    
    public void setSliderVerticalOrientation(boolean b){
        if(b){getSlider().setOrientation(DecimalSlider.VERTICAL);}
    }
    
    public boolean getEditValues(){ return editValues;}
    public void setEditValues(boolean b){
       editValues = b;
       textField.setEditable(editValues);
    }
    
    public int getPrecision(){return slider.getPrecision();}
    public void setPrecision(int n) {
        if(n>5){column = n;}
        slider.setPrecision(n);
    }    
    
    public double getValue(){return slider.getDecimalSliderValue();}    
    public void setValue(double d){
        slider.setDecimalSliderValue(d);
        textField.setText(String.valueOf(d));
    }
    
    /**
     * Returns the GUI element for display in the simulation.
     */
    public java.awt.Component graphic(Object obj) {
//        if(showValues){ return panel;
//        } else {return slider;}
        return panel;
    }
    
    
    /**
     * Sets the value of a descriptive label using the meter's label and the unit's symbol (abbreviation).
     */
    private void setLabelDefault() {
        String suffix = (unit.symbol().length() > 0) ? " ("+unit.symbol()+")" : "";
        if(modifyAction != null) 
            setLabel(StringUtility.capitalize(modifyAction.getLabel())+suffix);
    }

    /**
     * Sets the value of a descriptive label using the given string.
     */
    public void setLabel(String s){
        label = s;
        if(s == null || s.equals("") || !showBorder) panel.setBorder(new javax.swing.border.EmptyBorder(0,0,0,0));
        else panel.setBorder(new javax.swing.border.TitledBorder(s));
    }
    
    /**
     * @return the current instance of the descriptive label.
     */
    public String getLabel() {return label;}
    
    public void setShowBorder(boolean b) {
        showBorder = b;
        setLabel(label);
    }
    public boolean isShowBorder() {return showBorder;}
    
    /**
     * @return a handle to the DecimalSlider instance used by this slider device
     */
    public DecimalSlider getSlider() {return slider;}
    
    /**
     * @return a handle to the JTextField instance used by this slider device
     */    
    public JTextField getTextField() {return textField;}
    
    /**
     * @return a handle to the JPanel instance used by this slider device
     */    
    public JPanel getPanel() { return panel; }    
    /**
     * Performs no action (slider can not be accessed for setting).
     * Method exists so that slider and its properties (which can be set) will
     * show up on the property sheet.
     */
     public void setSlider(JSlider s) {}
     
    /**
     * Method to demonstrate and test the use of this class.  
     * Slider is used to control the diameter of a hard-sphere in MD simulation
     */
//    public static void main(String[] args) {
//         
//// inner class of Modifier that works with Slider*************************************************************
//           class DiameterModifier extends etomica.ModifierAbstract {
//
//                public double fullDiameter;
//                public double currentValue;
//                public etomica.graphics.DisplayPhase display;
//                private etomica.P2SquareWell p2SquareWell;
//                private etomica.SpeciesSpheresMono speciesSpheres01;
//                private JPanel valuePanel;
//
//                public DiameterModifier(etomica.P2SquareWell pot, etomica.SpeciesSpheresMono ssm){
//                    p2SquareWell = pot;
//                    speciesSpheres01 = ssm;
//                    fullDiameter = pot.getCoreDiameter();
//                }
//                public etomica.units.Dimension getDimension() {
//                    return etomica.units.Dimension.LENGTH; 
//                }
//                public double getValue() {
//                    return p2SquareWell.getCoreDiameter();
//                }
//                public void setValue(double d) { 
//                    p2SquareWell.setCoreDiameter(d);
//                    speciesSpheres01.setDiameter(d);
//                    display.repaint();
//                }
//                public void setDisplay(etomica.graphics.DisplayPhase display){
//                    this.display = display;
//                }
//            } 
////end of DiameterModifier class*********************************************************************
//           
//        Simulation.instance = new etomica.graphics.SimulationGraphic(); 
//
//        Phase phase0  = new Phase();
//            phase0.setLrcEnabled(false);
//        IntegratorHard integrator = new IntegratorHard();    
//        P2SquareWell p2SquareWell  = new P2SquareWell();
//        Controller controller0  = new Controller();
//        SpeciesSpheresMono speciesSpheres0  = new SpeciesSpheresMono();
//        DisplayPhase displayPhase0  = new DisplayPhase();
//        DeviceTrioControllerButton button = new DeviceTrioControllerButton();
////           button.setTrioControllerButtonShape("VERTICAL");     
//
//        DiameterModifier diaModifier= new DiameterModifier(p2SquareWell, speciesSpheres0);
//           diaModifier.setDisplay(displayPhase0);
//           
//        DeviceSlider mySlider = new DeviceSlider();
////           mySlider.setShowMinorValues(true);
//           mySlider.setPrecision(3);  // default "0" - working with integer without setting precision, better higher precesion than minimum and maximum
//           mySlider.setMinimum(0);  // possible to give double value
//           mySlider.setMaximum(6.2); // possible to give double value
//           mySlider.setModifier(diaModifier); // call modifier instance after setting precision, minimum, and maximum
////           mySlider.setValue(2*Default.atomSize);// if modifier instance is called first then setValue method should be called to set slider value
////           mySlider.setLabel("  ");  // without this, shows modifier demension  
//           mySlider.setShowValues(true); // default "false" - true makes panel to put Slider and TextField, which shows the values of slider
//           mySlider.setEditValues(true); // defaulst " false" - decide to edit the values after true setShowValues 
//           mySlider.setLabel("Diameter");
///*           mySlider.getJPanel().setBorder(new javax.swing.border.TitledBorder( // default border is null
//                                            new javax.swing.border.EtchedBorder(
//                                                javax.swing.border.EtchedBorder.RAISED, java.awt.Color.black, java.awt.Color.gray) 
//                                                ,""
//                                                ,javax.swing.border.TitledBorder.LEFT
//                                                ,javax.swing.border.TitledBorder.TOP
//                                                ,new java.awt.Font(null,java.awt.Font.BOLD,15)
//                                                ,java.awt.Color.black));   */         
////            mySlider.setSliderValueShape("HORIZONTAL"); // default "VERTICAL" 
////            mySlider.setSliderVerticalOrientation(true); // true is for vertically standing slider 
////            mySlider.getTextField().setHorizontalAlignment(mySlider.getTextField().LEFT); // default "CENTER"           
//                    
//        Simulation.instance.elementCoordinator.go();
//        SimulationGraphic.makeAndDisplayFrame(Simulation.instance);
//
//    }
}
