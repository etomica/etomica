package etomica.normalmode;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import etomica.api.IAction;
import etomica.data.DataInfo;
import etomica.data.IDataSink;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTable;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayTable;
import etomica.graphics.SimulationGraphic;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.units.Null;

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
		
		
		/*
		 * N atom Slider
		 */
		final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setBox(sim.box);
        nSlider.setSpecies(sim.species);
        nSlider.setMinimum(2);
        nSlider.setMaximum(40);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        
        int m = sim.nm.getOmegaSquared(sim.box).length;
        double[] omega2 = new double[m];
        for (int i=0; i<m; i++){
        	omega2[i] = sim.nm.getOmegaSquared(sim.box)[i][0];
        }
        
        data = new DataDoubleArray[1];
        data[0] = new DataDoubleArray(new int[]{m},omega2);
        
        omega2Table = new DataTable(data);
        
        nSlider.setPostAction(new IAction() {
        	
       	  	public void actionPerformed() {
       	  		int n = (int)nSlider.getValue(); 
       	  		System.out.println("n: "+n);     	  		
                if(n == 0) {
                	sim.integrator.setThermostatInterval(100);
                }
                else {
                	sim.integrator.setThermostatInterval(100/n);
                }
                
                if (oldN != n) {
                	Boundary boundary = new BoundaryRectangularPeriodic(sim.getSpace(), sim.getRandom(), n/sim.density);
                	sim.box.setBoundary(boundary);
                	sim.waveVectorFactory.makeWaveVectors(sim.box);
                	sim.coordinateDefinition.initializeCoordinates(new int[]{(int)nSlider.getValue()});
                }            
                
                oldN = n;
                try {
                	sim.integrator.reset();
     
                	sim.integrator.setWaveVectors(sim.waveVectorFactory.getWaveVectors());
                    sim.integrator.setWaveVectorCoefficients(sim.waveVectorFactory.getCoefficients());
                    sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(sim.box), sim.waveVectorFactory.getCoefficients());
                    sim.integrator.setEigenVectors(sim.nm.getEigenvectors(sim.box));
                    
                    int m = sim.nm.getOmegaSquared(sim.box).length;
                    double[] omega2 = new double[m];
                    for (int i=0; i<m; i++){
                    	omega2[i] = sim.nm.getOmegaSquared(sim.box)[i][0];
                    	System.out.print(omega2[i] +" ");
                    }
                    System.out.println();
                    
                    
                    data[0] = new DataDoubleArray(new int[]{m},omega2);
                    omega2Table = new DataTable(data);
                                       
                    DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Omega2", Null.DIMENSION, new int[]{m});
                    DataInfo dataInfo = new DataInfoTable("Omega2", new DataInfoDoubleArray[]{columnInfo}, m, null);
                    sink.putDataInfo(dataInfo);
                    sink.putData(omega2Table);
                 
                }
                
                catch (ConfigurationOverlapException e) {
                    throw new RuntimeException(e);
                }   	    
                
                resetAction.actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
       	  	int oldN = sim.box.getMoleculeList().getMoleculeCount();
        });
        
        add(nSlider);
        //end of N Slider
 
        //
        displayTable = new DisplayTable();
        sink = displayTable.getDataTable().makeDataSink();
        
        displayTable.setTransposed(false);
	    //displayTable.setColumnHeader(new DataTag[]{},"Omega2");
        
        DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Omega2", Null.DIMENSION, new int[]{m});
        DataInfo dataInfo = new DataInfoTable("Omega2", new DataInfoDoubleArray[]{columnInfo}, m, null);
        sink.putDataInfo(dataInfo);
        sink.putData(omega2Table);
        
        getPanel().tabbedPane.add("Omega^2", displayTable.graphic());
        
        //
        
        
        
        /*
		 * Temperature Slider
		 */
		temperatureSetter = new DeviceThermoSlider(sim.getController());
		temperatureSetter.setPrecision(1);
		temperatureSetter.setMinimum(0.0);
		temperatureSetter.setMaximum(10.0);
		temperatureSetter.setSliderMajorValues(5);
		temperatureSetter.setIntegrator(sim.integrator);
		temperatureSetter.setIsothermal();
		temperatureSetter.setTemperature(sim.temperature);
		
		
		final IAction temperatureAction = new IAction() {
            public void actionPerformed() {
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
		Space sp = Space.getInstance(1);
		NormalModeAnalysisDisplay1DGraphic simGraphic = new NormalModeAnalysisDisplay1DGraphic(new NormalModeAnalysisDisplay1D(sp), sp);
		SimulationGraphic.makeAndDisplayFrame(simGraphic.getPanel(), APP_NAME);
		
	}
	
	public static class Applet extends javax.swing.JApplet {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public void init(){
			getRootPane().putClientProperty(APP_NAME, Boolean.TRUE);
			Space sp = Space.getInstance(1);
			NormalModeAnalysisDisplay1DGraphic nm1Dgraphic = new NormalModeAnalysisDisplay1DGraphic(new NormalModeAnalysisDisplay1D(sp), sp);
			getContentPane().add(nm1Dgraphic.getPanel());
		}
	}
	
	private DeviceThermoSlider temperatureSetter; 
	private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "1-D Harmonic Oscillator";
	private static final int REPAINT_INTERVAL = 10;
	protected NormalModeAnalysisDisplay1D sim;
	protected final IAction resetAction;
	protected DataTable omega2Table;
	protected DataDoubleArray[] data;
	protected final DisplayTable displayTable;
	protected IDataSink sink;

}
