package etomica.normalmode;

import etomica.api.IAction;
import etomica.api.IData;
import etomica.data.DataInfo;
import etomica.data.IDataSink;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTable;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeRandom;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayTable;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.units.Energy;
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
		
		temperatureSetter.setSliderPostAction(new IAction() {
            public void actionPerformed() {
            	
            	sim.nm.setTemperature(sim.integrator.temperature);
                
                int m = sim.nm.getOmegaSquared(sim.box).length;
                double[] omega2 = new double[m];
                
                for (int i=0; i<m; i++){
                	omega2[i] = sim.nm.getOmegaSquared(sim.box)[i][0];
                	stringWV[i]=String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(0));
                }
                
                data[0] = new DataDoubleArray(new int[]{m},omega2);
                omega2Table = new DataTable(data);
               
                /*
                 * Harmonic Energy
                 */
                AHarm = new DataDouble();
                double AHarmonic =0;
                
                coeffs = sim.nm.getWaveVectorFactory().getCoefficients();
                for(int i=0; i<omega2.length; i++) {
                        if (!Double.isInfinite(omega2[i])) {
                            AHarmonic += coeffs[i] * Math.log(omega2[i]*coeffs[i] /
                                    (sim.integrator.temperature*Math.PI));
                    }
                }
                AHarm.E(AHarmonic);
                displayAHarmonic.putDataInfo(dataInfoA);
                displayAHarmonic.putData(AHarm);
                displayAHarmonic.repaint();
                
                DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Omega^2", Null.DIMENSION, new int[]{m});
                DataInfo dataInfo = new DataInfoTable("Omega^2", new DataInfoDoubleArray[]{columnInfo}, m, stringWV);
                sink.putDataInfo(dataInfo);
                sink.putData(omega2Table);
            	                
                resetAction.actionPerformed();
            	
		    }
		});
       
		
		// end of Temperature Slider
		
		
		
		
		
		
		
		/*
		 * N atom Slider
		 */
		final DeviceNSelector nSlider = new DeviceNSelector(sim.getController());
        nSlider.setBox(sim.box);
        nSlider.setSpecies(sim.species);
        nSlider.setMinimum(2);
        nSlider.setMaximum(50);
        nSlider.setLabel("Number of Atoms");
        nSlider.setShowBorder(true);
        nSlider.setShowValues(true);
        
        int m = sim.nm.getOmegaSquared(sim.box).length;
        double[] omega2 = new double[m];
        stringWV = new String[m];
     
        for (int i=0; i<m; i++){
        	omega2[i] = sim.nm.getOmegaSquared(sim.box)[i][0];
        	stringWV[i]=String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(0));
        }
        
        data = new DataDoubleArray[1];
        data[0] = new DataDoubleArray(new int[]{m},omega2);
        omega2Table = new DataTable(data);
        
        /*
         * Harmonic Energy
         */
        
        AHarm = new DataDouble();
        double AHarmonic =0;
        
        coeffs = sim.nm.getWaveVectorFactory().getCoefficients();
        
        for(int i=0; i<omega2.length; i++) {
                if (!Double.isInfinite(omega2[i])) {
                    AHarmonic += coeffs[i] * Math.log(omega2[i]*coeffs[i] /
                            (sim.integrator.temperature*Math.PI));
            }
        }
        AHarm.E(AHarmonic);
        
        nSlider.setPostAction(new IAction() {
        	
       	  	public void actionPerformed() {
       	  		int n = (int)nSlider.getValue();                 
       	  	
                if (oldN != n) {
                	Boundary boundary = new BoundaryRectangularPeriodic(sim.getSpace(), sim.getRandom(), n/sim.density);
                	sim.box.setBoundary(boundary);
                	sim.waveVectorFactory.makeWaveVectors(sim.box);
                	sim.coordinateDefinition.initializeCoordinates(new int[]{(int)nSlider.getValue()});
                }            
                
                oldN = n;
                try {
                	sim.integrator.reset();
                	
                	sim.nm.setTemperature(sim.integrator.temperature);
                	sim.integrator.setWaveVectors(sim.waveVectorFactory.getWaveVectors());
                    sim.integrator.setWaveVectorCoefficients(sim.waveVectorFactory.getCoefficients());
                    sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(sim.box), sim.waveVectorFactory.getCoefficients());
                    sim.integrator.setEigenVectors(sim.nm.getEigenvectors(sim.box));
                    
                    int m = sim.nm.getOmegaSquared(sim.box).length;
                    double[] omega2 = new double[m];
                    
                    for (int i=0; i<m; i++){
                    	omega2[i] = sim.nm.getOmegaSquared(sim.box)[i][0];
                    	stringWV[i]=String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(0));
                    }
                    
                    data[0] = new DataDoubleArray(new int[]{m},omega2);
                    omega2Table = new DataTable(data);
                   
                    
                    /*
                     * Harmonic Energy
                     */
                    AHarm = new DataDouble();
                    
                    double AHarmonic =0;
                    coeffs = sim.nm.getWaveVectorFactory().getCoefficients();
                    
                    for(int i=0; i<omega2.length; i++) {
                            if (!Double.isInfinite(omega2[i])) {
                                AHarmonic += coeffs[i] * Math.log(omega2[i]*coeffs[i] /
                                        (sim.integrator.temperature*Math.PI));
                        }
                    }
                    AHarm.E(AHarmonic);
                    DataInfoDouble dataInfoA = new DataInfoDouble("AHarmonic", Energy.DIMENSION);
                    displayAHarmonic.putDataInfo(dataInfoA);
                    displayAHarmonic.putData(AHarm);
                    displayAHarmonic.repaint();
                    
                    DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Omega^2", Null.DIMENSION, new int[]{m});
                    DataInfo dataInfo = new DataInfoTable("Omega^2", new DataInfoDoubleArray[]{columnInfo}, m, stringWV);
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
        
        //end of N Slider
        

        displayAHarmonic = new DisplayTextBox();
        displayAHarmonic.setPrecision(10);
        dataInfoA = new DataInfoDouble("AHarmonic", Energy.DIMENSION);
        displayAHarmonic.putDataInfo(dataInfoA);
        displayAHarmonic.putData(AHarm);
        displayAHarmonic.setLabel("Harmonic Energy");
      
        add(nSlider);
        add(temperatureSetter);
    
		
        
        getDisplayBox(sim.box).setColorScheme(new ColorSchemeRandom(sim, sim.box,sim.getRandom()));
			
        
        /*
         * tabbed-pane for wavevectors with corresponding omega2
         */
        displayTable = new DisplayTable();
        sink = displayTable.getDataTable().makeDataSink();
        
        displayTable.setTransposed(false);
        
        DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Omega^2", Null.DIMENSION, new int[]{m});
        DataInfo dataInfoTable = new DataInfoTable("Omega^2", new DataInfoDoubleArray[]{columnInfo}, m, stringWV);
        sink.putDataInfo(dataInfoTable);
        sink.putData(omega2Table);
        
        getPanel().tabbedPane.add("Omega^2", displayTable.graphic());
        
        //
        add(displayAHarmonic);
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
	protected String[] stringWV;
	protected double[] coeffs;
	protected IData AHarm;
	protected DisplayTextBox displayAHarmonic;
	protected DataInfoDouble dataInfoA; 

}
