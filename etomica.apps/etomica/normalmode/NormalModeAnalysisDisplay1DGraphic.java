package etomica.normalmode;

import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import etomica.api.IAction;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.units.Temperature;

/**
 * 
 * 
 * @author Tai Boon Tan
 */
public class NormalModeAnalysisDisplay1DGraphic extends SimulationGraphic {


	public NormalModeAnalysisDisplay1DGraphic(final NormalModeAnalysisDisplay1D simulation, Space space) {
		
		super(simulation, TABBED_PANE, APP_NAME,REPAINT_INTERVAL, space, simulation.getController());
		this.sim = simulation;
		
		resetAction = getController().getSimRestart().getDataResetAction();
		
		final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setSpecies(sim.species);
        nSlider.setBox(sim.box);
        nSlider.setMinimum(0);
        nSlider.setMaximum(40);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        
        nSlider.setPostAction(new IAction() {
            public void actionPerformed() {
                int n = (int)nSlider.getValue();
                if(n == 0) {
                	sim.integrator.setThermostatInterval(50);
                }
                else {
                	sim.integrator.setThermostatInterval(50/n);
                }
                
                if (oldN < n) {
                	Boundary boundary = new BoundaryRectangularPeriodic(sim.getSpace(), sim.getRandom(), n/sim.density);
                	System.out.println("numAtom: "+ n);
                	sim.box.setBoundary(boundary);
                	CoordinateDefinition coordinateDefinition = new CoordinateDefinitionLeaf(sim, sim.box, sim.primitive, sim.getSpace());
                    coordinateDefinition.initializeCoordinates(new int[]{(int)nSlider.getValue()});
                }
                oldN = n;
                try {
                    sim.integrator.reset();
                }
                catch (ConfigurationOverlapException e) {
                    throw new RuntimeException(e);
                }
                resetAction.actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
                      
            int oldN = sim.box.getMoleculeList().getAtomCount();
           
        });
        
        add(nSlider);
       
		/*
		 * Temperature Slider
		 */
		temperatureSetter = new DeviceThermoSlider(sim.getController());
		temperatureSetter.setPrecision(1);
		temperatureSetter.setIsothermal();
		temperatureSetter.setMinimum(0.0);
		temperatureSetter.setMaximum(5.0);
		temperatureSetter.setTemperature(0.0);
		temperatureSetter.setSliderMajorValues(4);
		temperatureSetter.setUnit(Temperature.SIM_UNIT);
		temperatureSetter.setIntegrator(sim.integrator);
		
        final IAction temperatureAction = new IAction() {
            public void actionPerformed() {
                resetAction.actionPerformed();
		    }
		};
		ActionListener isothermalListener = new ActionListener() {
		    public void actionPerformed(ActionEvent event) {
		        temperatureAction.actionPerformed();
		    }
		};

		temperatureSetter.setSliderPostAction(temperatureAction);
        temperatureSetter.addRadioGroupActionListener(isothermalListener);
		
		add(temperatureSetter);
		// end of Temperature Slider
	}
	
	public static void main(String[] args){
		
		NormalModeAnalysisDisplay1D sim = new NormalModeAnalysisDisplay1D();
		NormalModeAnalysisDisplay1DGraphic simGraphic = new NormalModeAnalysisDisplay1DGraphic(sim, Space.getInstance(1));
		SimulationGraphic.makeAndDisplayFrame(simGraphic.getPanel(), APP_NAME);
		
	}
	


	
	private DeviceThermoSlider temperatureSetter; 
	private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "1-D Harmonic Oscillator";
	private static final int REPAINT_INTERVAL = 10;
	protected NormalModeAnalysisDisplay1D sim;
	protected final IAction resetAction;

}
