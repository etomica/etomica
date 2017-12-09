/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;
import etomica.data.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataTable;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.listener.IntegratorListenerAction;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Null;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * 
 * 
 * @author Tai Boon Tan
 */
public class NormalModeAnalysisDisplay2DGraphic extends SimulationGraphic {

	public NormalModeAnalysisDisplay2DGraphic(final NormalModeAnalysisDisplay2D simulation, Space space) {
		
		super(simulation, TABBED_PANE, APP_NAME,REPAINT_INTERVAL);
		this.sim = simulation;
		
		DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
		
		
		/*
		 * harmonic energy                                                                            
		 */
		MeterHarmonicEnergy heMeter = new MeterHarmonicEnergy(sim.coordinateDefinition, sim.nm);
		
		AccumulatorHistory heHistory = new AccumulatorHistory();
		heHistory.setTimeDataSource(timeCounter);
	    
		final AccumulatorAverageCollapsing heAccumulator = new AccumulatorAverageCollapsing();
        heAccumulator.setPushInterval(10);
        DataFork heFork = new DataFork(new IDataSink[]{heHistory, heAccumulator});
        DataPump hePump = new DataPump(heMeter, heFork);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(hePump);
        pumpListener.setInterval(60);
        sim.integrator.getEventManager().addListener(pumpListener);
        heHistory.setPushInterval(5);
		
        DisplayPlot ePlot = new DisplayPlot();
        heHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.setLegend(new DataTag[]{heHistory.getTag()}, "Harmonic Energy");
        
        ePlot.getPlot().setTitle("Energy History");
        ePlot.setDoLegend(true);
        ePlot.setLabel("Energy");
        
        final DisplayTextBoxesCAE heDisplay = new DisplayTextBoxesCAE();
        heDisplay.setAccumulator(heAccumulator);
   
        
        /*
		 * Temperature Slider
		 */
		temperatureSetter = new DeviceThermoSlider(sim.getController(), sim.integrator);
		temperatureSetter.setIsothermalButtonsVisibility(false);
		temperatureSetter.setPrecision(1);
		temperatureSetter.setMinimum(0.0);
		temperatureSetter.setMaximum(10.0);
		temperatureSetter.setSliderMajorValues(5);
		temperatureSetter.setIsothermal();
		temperatureSetter.setTemperature(sim.temperature);
		
		temperatureSetter.setSliderPostAction(new IAction() {
            public void actionPerformed() {
            	
            	int m = sim.nm.getOmegaSquared().length;
                double[] omega2 = new double[m];
                
                waveVectorSlider.setMaximum(m);
                waveVectorSlider.setIntegrator(sim.integrator);
                //Array
                if(sim.integrator.getWaveVectorNum() >= m){
                	waveVectorSlider.setMaximum(m);
                	sim.integrator.setWaveVectorNum(m-1);
                }
                
                
                if (sim.integrator.isOneWV()){
                	int wvNumUsed = sim.integrator.getWaveVectorNum();
                	for (int i=0; i<m; i++){
                    	omega2[i] = sim.nm.getOmegaSquared()[i][0];
                    	
                    	if(i==wvNumUsed){
                    		stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                        	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";
                    	} else {
                    		stringWV[i]= String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                        	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1));
                    	}
                    }
                } else {
                
                    for (int i=0; i<m; i++){
                    	omega2[i] = sim.nm.getOmegaSquared()[i][0];
                    	stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                    	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";
                    }
                }
                data[0] = new DataDoubleArray(new int[]{m},omega2);
                omega2Table = new DataTable(data);
               
                /*
                 * Harmonic Free Energy
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
            	                
                getController().getSimRestart().getDataResetAction().actionPerformed();
            	
		    }
		});
       
		
		// end of Temperature Slider
		
	
		
		
		/*
		 * N Cell Slider
		 */
		final DeviceCellNumXYSlider cellSlider = new DeviceCellNumXYSlider(sim.getController());
		cellSlider.setBox(sim.box);
		cellSlider.setSpecies(sim.species);
		cellSlider.setMaximum(sim.getDimx());
        
        int m = sim.nm.getOmegaSquared().length;
        double[] omega2 = new double[m];
        stringWV = new String[m];
     
        if (sim.integrator.isOneWV()){
        	int wvNumUsed = sim.integrator.getWaveVectorNum();
        	for (int i=0; i<m; i++){
            	omega2[i] = sim.nm.getOmegaSquared()[i][0];
            	
            	if(i==wvNumUsed){
            		stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";
            	} else {
            		stringWV[i]= String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1));
            	}
            }
        } else {
        
            for (int i=0; i<m; i++){
            	omega2[i] = sim.nm.getOmegaSquared()[i][0];
            	stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
            	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";
            }
        }
        
        data = new DataDoubleArray[1];
        data[0] = new DataDoubleArray(new int[]{m},omega2);
        omega2Table = new DataTable(data);
        
        /*
         * Harmonic Free Energy
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
        
        cellSlider.setXSliderPostAction(new IAction() {
        	
       	  	public void actionPerformed() {
       	  		
       	  		int nx = (int)cellSlider.getXCellNum(); 
       	  	
                if (oldNx != nx ) {
                	int[] nCells = new int[]{nx, (int)cellSlider.getYCellNum()};
                	
                	Vector[] dimension = sim.space.makeVectorArray(sim.space.D());
                    for (int i=0; i<sim.space.D(); i++){
                    	dimension[i].Ea1Tv1(nCells[i], sim.primitive.vectors()[i]);
                    }
                	                	       	
                	Boundary boundary = new BoundaryDeformablePeriodic(sim.getSpace(), dimension);
                	sim.box.setBoundary(boundary);
                	
                	sim.waveVectorFactory.makeWaveVectors(sim.box);
                	sim.coordinateDefinition.initializeCoordinates(nCells);
                	
                }            
                
                oldNx = nx;

            	sim.integrator.reset();
            	
            	sim.integrator.setWaveVectors(sim.waveVectorFactory.getWaveVectors());
                sim.integrator.setWaveVectorCoefficients(sim.waveVectorFactory.getCoefficients());
                sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(), sim.waveVectorFactory.getCoefficients());
                sim.integrator.setEigenVectors(sim.nm.getEigenvectors());
                                                
                int m = sim.nm.getOmegaSquared().length;
                double[] omega2 = new double[m];
                
                waveVectorSlider.setMaximum(m);
                waveVectorSlider.setIntegrator(sim.integrator);
                //Array
                if(sim.integrator.getWaveVectorNum() >= m){
                	waveVectorSlider.setMaximum(m);
                	sim.integrator.setWaveVectorNum(m-1);
                }
                
                
                if (sim.integrator.isOneWV()){
                	int wvNumUsed = sim.integrator.getWaveVectorNum();
                	for (int i=0; i<m; i++){
                    	omega2[i] = sim.nm.getOmegaSquared()[i][0];
                    	
                    	if(i==wvNumUsed){
                    		stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                        	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";
                    	} else {
                    		stringWV[i]= String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                        	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1));
                    	}
                    }
                } else {
                
                    for (int i=0; i<m; i++){
                    	omega2[i] = sim.nm.getOmegaSquared()[i][0];
                    	stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                    	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";
                    }
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

                getController().getSimRestart().getDataResetAction().actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
       	  	int oldNx = (int)cellSlider.getXCellNum();
        });
        
        //end of XCell Post-Action
        
        cellSlider.setYSliderPostAction(new IAction() {
        	
       	  	public void actionPerformed() {
       	  		
       	  		int ny = (int)cellSlider.getYCellNum(); 
       	  	
                if (oldNy != ny ) {
                	int[] nCells = new int[]{(int)cellSlider.getXCellNum(), ny};
                	
                	Vector[] dimension = sim.space.makeVectorArray(sim.space.D());
                    for (int i=0; i<sim.space.D(); i++){
                    	dimension[i].Ea1Tv1(nCells[i], sim.primitive.vectors()[i]);
                    }
                	                	       	
                	Boundary boundary = new BoundaryDeformablePeriodic(sim.getSpace(), dimension);
                	sim.box.setBoundary(boundary);
                	sim.waveVectorFactory.makeWaveVectors(sim.box);
                	sim.coordinateDefinition.initializeCoordinates(nCells);
                	
                }            
                
                oldNy = ny;
                try {
                	sim.integrator.reset();
                	
                	sim.integrator.setWaveVectors(sim.waveVectorFactory.getWaveVectors());
                    sim.integrator.setWaveVectorCoefficients(sim.waveVectorFactory.getCoefficients());
                    sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(), sim.waveVectorFactory.getCoefficients());
                    sim.integrator.setEigenVectors(sim.nm.getEigenvectors());
                                                    
                    int m = sim.nm.getOmegaSquared().length;
                    double[] omega2 = new double[m];
                    
                    waveVectorSlider.setMaximum(m);
                    waveVectorSlider.setIntegrator(sim.integrator);
                    //Array
                    if(sim.integrator.getWaveVectorNum() >= m){
                    	waveVectorSlider.setMaximum(m);
                    	sim.integrator.setWaveVectorNum(m-1);
                    }
                    
                    
                    if (sim.integrator.isOneWV()){
                    	int wvNumUsed = sim.integrator.getWaveVectorNum();
                    	for (int i=0; i<m; i++){
                        	omega2[i] = sim.nm.getOmegaSquared()[i][0];
                        	
                        	if(i==wvNumUsed){
                        		stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                            	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";
                        	} else {
                        		stringWV[i]= String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                            	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1));
                        	}
                        }
                    } else {
                    
	                    for (int i=0; i<m; i++){
	                    	omega2[i] = sim.nm.getOmegaSquared()[i][0];
	                    	stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
	                    	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";
	                    }
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
                
                getController().getSimRestart().getDataResetAction().actionPerformed();
                getDisplayBox(sim.box).repaint();
            }
       	 int oldNy = (int)cellSlider.getYCellNum();
        });
        
        //end of NCell Slider
        
        
        /*
         *  Wave vectors Slider
         */ 
           
        waveVectorSlider = new DeviceWaveVectorSlider(sim.getController());
        waveVectorSlider.setMinimum(0);
        waveVectorSlider.setMaximum(m);
        waveVectorSlider.setOneWVButtonsVisibility(false);
        waveVectorSlider.setIntegrator(sim.integrator);
        
        final IAction waveVectorAction = new IAction() {
        	
        	public void actionPerformed() {
        		
                int m = sim.nm.getOmegaSquared().length;
                double[] omega2 = new double[m];
                
                int wvNumUsed = (int)waveVectorSlider.getWaveVectorNum();
                
                if (waveVectorSlider.isOneWV()){
	                for (int i=0; i<m; i++){
	                	
	                	omega2[i] = sim.nm.getOmegaSquared()[i][0];
	                	if (i==wvNumUsed){
	                		stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
	                    	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";
	                		
	                	} else {
	                	
	                		stringWV[i]=String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
	                    	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1));
	                	}
	                }
                } else {
                	for (int i=0; i<m; i++){
                    	omega2[i] = sim.nm.getOmegaSquared()[i][0];
                    	stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(0))+", "
                    	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].getX(1))+">";	
                    }
                }
                
                
                data[0] = new DataDoubleArray(new int[]{m},omega2);
                omega2Table = new DataTable(data);
                
                DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Omega^2", Null.DIMENSION, new int[]{m});
                DataInfo dataInfo = new DataInfoTable("Omega^2", new DataInfoDoubleArray[]{columnInfo}, m, stringWV);
                sink.putDataInfo(dataInfo);
                sink.putData(omega2Table);
            	                
                getController().getSimRestart().getDataResetAction().actionPerformed();
                
		    }
		};
		
		ActionListener isOneWVListener = new ActionListener(){
			public void actionPerformed (ActionEvent event){
				waveVectorAction.actionPerformed();
				
			}
		};
        
		waveVectorSlider.setSliderPostAction(waveVectorAction);
		waveVectorSlider.addRadioGroupActionListener(isOneWVListener);
        // end wave vectors slider
        
        
        displayAHarmonic = new DisplayTextBox();
        displayAHarmonic.setPrecision(10);
        dataInfoA = new DataInfoDouble("AHarmonic", Energy.DIMENSION);
        displayAHarmonic.putDataInfo(dataInfoA);
        displayAHarmonic.putData(AHarm);
        displayAHarmonic.setLabel("Harmonic Free Energy");
      
  
        /*
         * tabbed-pane for wavevectors with corresponding omega2
         */
        displayTable = new DisplayTable();

        displayTable.setTransposed(false);
        
        DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Omega^2", Null.DIMENSION, new int[]{m});
        DataInfo dataInfoTable = new DataInfoTable("Omega^2", new DataInfoDoubleArray[]{columnInfo}, m, stringWV);
        sink = displayTable.getDataTable().makeDataSink(dataInfoTable);
        sink.putDataInfo(dataInfoTable);
        sink.putData(omega2Table);
        
        getPanel().tabbedPane.add("Omega^2", displayTable.graphic());
        
        //
        
        
        getController().getDataStreamPumps().add(hePump);
        
        IAction resetAction = new IAction(){
        	public void actionPerformed(){
        		heDisplay.putData(heAccumulator.getData());
        		heDisplay.repaint();
        		
        		getDisplayBox(sim.box).graphic().repaint();
        	}
        };
        
        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);
        
        add(cellSlider);
        add(temperatureSetter);
        add(waveVectorSlider);
        add(displayAHarmonic);
        add(ePlot);
        add(heDisplay);
	}
	
	
	
	public static void main(String[] args){
		Space sp = Space.getInstance(2);
		NormalModeAnalysisDisplay2DGraphic simGraphic = new NormalModeAnalysisDisplay2DGraphic(new NormalModeAnalysisDisplay2D(sp), sp);
		SimulationGraphic.makeAndDisplayFrame(simGraphic.getPanel(), APP_NAME);
		
	}
	
	public static class Applet extends javax.swing.JApplet {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public void init(){
			getRootPane().putClientProperty(APP_NAME, Boolean.TRUE);
			Space sp = Space.getInstance(2);
			NormalModeAnalysisDisplay2DGraphic nm2Dgraphic = new NormalModeAnalysisDisplay2DGraphic(new NormalModeAnalysisDisplay2D(sp), sp);
			getContentPane().add(nm2Dgraphic.getPanel());
		}
	}
	
	private DeviceThermoSlider temperatureSetter; 
	private DeviceWaveVectorSlider waveVectorSlider;
	private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "2-D Harmonic Oscillator";
	private static final int REPAINT_INTERVAL = 10;
	protected NormalModeAnalysisDisplay2D sim;
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
