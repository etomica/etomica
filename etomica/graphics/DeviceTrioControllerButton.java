package etomica.graphics;
import etomica.*;
import etomica.action.*;
import javax.swing.JPanel;
import java.awt.event.*;

/**
 * Device comprising three buttons: (1) attaches to a controller to toggle its pause/resume state; 
 * (2)causes a reset (restart) action of simulation to be performed; (3) causes reset of meter averages
 * only.
 *
 * @author Sang Kyu Kwak
 */
 
 /* History of changes.
  * 07/03/02 Created
  * 01/04/03 Added methods to access buttons
  */
public class DeviceTrioControllerButton extends Device {
    
    private JPanel jp;
    private DeviceControllerButton button1;
    private DeviceButton button2;
    private DeviceButton button3; 
    private SimulationRestart simReset;
    private MeterReset meterAverageReset;
	private Simulation simulation;	
	private double width;
	private boolean firstResized = true;
	private String shape;

    
    public DeviceTrioControllerButton(){
        this(Simulation.instance);
    }
    
    /**
     * Constructor if instance of controller is not added.
     */
    public DeviceTrioControllerButton(Simulation sim) {
        super(sim);
        simulation = sim;
        simReset = new SimulationRestart(sim);
        meterAverageReset = new MeterReset(sim);
        
        jp = new JPanel(new java.awt.GridLayout(1, 3)); //default shape of panel
        jp.setBorder(new javax.swing.border.TitledBorder("Control"));
        jp.setOpaque(false);
//        jp.setBackground(DefaultGraphic.BORDER_COLOR);
/*        jp.setBorder(new javax.swing.border.TitledBorder(
                         new javax.swing.border.EtchedBorder(
                             javax.swing.border.EtchedBorder.RAISED, java.awt.Color.red, java.awt.Color.blue) 
                             ,"Control"
                             ,javax.swing.border.TitledBorder.LEFT
                             ,javax.swing.border.TitledBorder.TOP
                             ,new java.awt.Font(null,java.awt.Font.BOLD,15)
                             ,java.awt.Color.black));
                             */
        button1 = new DeviceControllerButton(sim.space);
        button2 = new DeviceButton(sim.space);
        button2.setAction(new ActionGraphic(simReset));
        
        button3 = new DeviceButton(sim.space);
        button3.setAction(new ActionGraphic(meterAverageReset));
        
        jp.add(button1.graphic()); 
        jp.add(button2.graphic());  
        jp.add(button3.graphic());
                
		
        class MyMediator extends Mediator.ControllerNull {
            MyMediator(Mediator m) {
                super(m);
                setSuperceding(true);
            }
            /**
            * Causes addition of button that toggles controller state.
            */
            public void add(Controller c) {
                if(DeviceTrioControllerButton.this.getController() == null) {
                    setController(c);
                }            
            }
        }//end of MyMediator
        
        sim.elementCoordinator.addMediatorPair(new MyMediator(sim.elementCoordinator)); 
        
        setShape("VERTICAL");
    }
    
    /**
     * Constructor if instance of controller is added.
     * Inner class ControllerNull.NoAction causes previously added controller button to be ignored
     */        
    public DeviceTrioControllerButton(Simulation sim, Controller c) {
        this(sim);
        setController(c);        
        sim.elementCoordinator.addMediatorPair(new Mediator.ControllerNull.NoAction(sim.elementCoordinator));
    }
    
    // final due to being called in contructor
    /**
     * Sets the controller that is toggled by this device.
     */
    public final void setController(Controller c) {
        button1.setController(c);
    }
    
    /**
     * Returns the controller instance toggled by this device.
     */
    public Controller getController() {
        return button1.getController();
    }

    /**
     * Returns the JPanel that contains the buttons.
     */
    public java.awt.Component graphic(Object obj) {
        return jp;
    }
    
    public DeviceControllerButton controllerButton() {return button1;}
    public DeviceButton resetButton() {return button2;}
    public DeviceButton resetAveragesButton() {return button3;}
 
    /**
     * Sets controller toggle button to read "Start"
     */
    public void reset() {button1.reset();}
    
    /**
     * Sets display shape of the device.  The argument must be a string having
     * one of the following values:<ul>
     * <li>"HORIZONTAL" causes buttons to be laid out in a row.
     * <li>"VERTICAL" causes buttons to be laid out in a column (default)
     * <li>"AUTOMATIC" causes buttons to adopt layout based on size of available space
     * </ul>
     * Any other value is ignored.
     */
    public void setShape(String s){
        shape = s;
        if(s=="HORIZONTAL"){jp.setLayout(new java.awt.GridLayout(1,3));jp.updateUI();}
        if(s=="VERTICAL"){jp.setLayout(new java.awt.GridLayout(3,1));jp.updateUI();}
        if(s=="AUTOMATIC"){((SimulationGraphic)simulation).panel().addComponentListener(new ComponentEventControllerButton());}
    }
    public String getShape() {return shape;}
    
    /**
     * Inner class that catches action of simulation panel 
     */        
     private class ComponentEventControllerButton implements ComponentListener {
  
        public void componentHidden(ComponentEvent e){}
        public void componentMoved(ComponentEvent e){}
        public void componentShown(ComponentEvent e){}
        public void componentResized(ComponentEvent e){
            if(firstResized) { width = ((SimulationGraphic)simulation).panel().getWidth();firstResized =false;}
            if((double)((SimulationGraphic)simulation).panel().getWidth()<width){
                jp.setLayout(new java.awt.GridLayout(3,1));
                jp.updateUI();
            } else { 
                jp.setLayout(new java.awt.GridLayout(1,3)); 
                jp.updateUI();
            }
        }//end componentResized
    }//end ComponentEventControllerButton
 
    /**
     * main method to show how to work with this class 
     */        
     public static void main(String[] args) {
        
        Simulation.instance = new etomica.graphics.SimulationGraphic(); 

        Phase phase0  = new Phase();
            phase0.setLrcEnabled(false);
        IntegratorHard integrator = new IntegratorHard();    
        P2SquareWell p2Squarewell0  = new P2SquareWell();
        Controller controller0  = new Controller();
        SpeciesSpheresMono speciesSpheres0  = new SpeciesSpheresMono();
        DisplayPhase displayPhase0  = new DisplayPhase();
        MeterTemperature meterEnergy = new MeterTemperature();
        DisplayTable table = new DisplayTable();
            table.setWhichValues(new MeterAbstract.ValueType[] {MeterAbstract.CURRENT, MeterAbstract.AVERAGE});

        DeviceTrioControllerButton button = new DeviceTrioControllerButton();
            button.setShape("VERTICAL"); //three choices "HORIZONTAL", "AUTOMATIC"          
//        DeviceTrioControllerButton button = new DeviceTrioControllerButton(Simulation.instance, Simulation.instance.controller(0)); 
//          button.setShape("VERTICAL"); //three choices "HORIZONTAL", "AUTOMATIC"
        
        Simulation.instance.elementCoordinator.go();
        SimulationGraphic.makeAndDisplayFrame(Simulation.instance);
    }
    
       
}//end DeviceTrioControllerButton