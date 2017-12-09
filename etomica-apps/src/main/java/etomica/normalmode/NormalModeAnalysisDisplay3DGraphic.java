/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;
import etomica.data.*;
import etomica.data.history.HistoryScrolling;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataTable;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.graphics.*;
import etomica.listener.IntegratorListenerAction;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Null;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Harmonic Oscillator 3D
 * 
 * @author Tai Boon Tan
 */
public class NormalModeAnalysisDisplay3DGraphic extends SimulationGraphic {

	public NormalModeAnalysisDisplay3DGraphic(final NormalModeAnalysisDisplay3D simulation, Space space) {
		
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
        IntegratorListenerAction hePumpListener = new IntegratorListenerAction(hePump);
        hePumpListener.setInterval(60);
        sim.integrator.getEventManager().addListener(hePumpListener);
        heHistory.setPushInterval(5);

        final DisplayTextBoxesCAE heDisplay = new DisplayTextBoxesCAE();
        heDisplay.setAccumulator(heAccumulator);
        
        
        /*
         * harmonic energy as a function of atomic displacement Q
         */
        //HistoryScrolling is keeping track of the harmonic energy 
        final AccumulatorHistory heQHistory = new AccumulatorHistory(new HistoryScrolling(1)); 
        heFork.addDataSink(heQHistory);
   
        meterHarmonicCoordinate = new MeterHarmonicCoordinate(sim.coordinateDefinition);

        /*
         * soft-sphere potential energy as a function of atomic displacement Q
         */
		AccumulatorHistory peHistory = new AccumulatorHistory();
		peHistory.setTimeDataSource(timeCounter);
		
        final AccumulatorHistory peQHistory = new AccumulatorHistory(new HistoryScrolling(1)); 
        
        DataFork peFork = new DataFork(new IDataSink[]{peHistory, peQHistory});
        
        DataProcessor dataProcessor = new DataProcessor() {

            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
				dataInfo = new DataInfoDouble("Potential Energy",Energy.DIMENSION);
				data = new DataDouble();
								
				return dataInfo;
			}
		
			protected IData processData(IData inputData) {
				data.x = inputData.getValue(0)-sim.latticeEnergy;
				//System.out.println("Soft Sphere energy: "+ data.x);
				//System.exit(1);
				return data;
			}
			
			IDataInfo dataInfo;
			DataDouble data;
		};
		
        DataPump pePump = new DataPump(sim.meterPE, dataProcessor);
        dataProcessor.setDataSink(peFork);
        
        IntegratorListenerAction pePumpListener = new IntegratorListenerAction(pePump);
        pePumpListener.setInterval(60);
        sim.integrator.getEventManager().addListener(pePumpListener);
        peHistory.setPushInterval(5);
       
        /*
         * Potential Well Plot
         */
        final DisplayPlot eQPlot = new DisplayPlot();
        peQHistory.setDataSink(eQPlot.getDataSet().makeDataSink());
        heQHistory.setDataSink(eQPlot.getDataSet().makeDataSink());
        eQPlot.getPlot().setYLabel("Energy");
        eQPlot.setLabel("Potential Well");
        eQPlot.setDoDrawLines(new DataTag[]{peQHistory.getTag(), heQHistory.getTag()}, false);
        
        heQHistory.setTimeDataSource(meterHarmonicCoordinate);
        peQHistory.setTimeDataSource(meterHarmonicCoordinate);
        
        eQPlot.setDoClear(false);
        
		/*
		 * Energy Plot
		 */
        DisplayPlot ePlot = new DisplayPlot();
        peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        heHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        ePlot.getPlot().setYLabel("Energy");
        
        ePlot.getPlot().setTitle("Energy History");
        ePlot.setDoLegend(true);
        ePlot.setLabel("Energy History");
        
        /*
		 * Temperature Slider
		 */
		temperatureSetter = new DeviceThermoSlider(sim.getController(), sim.integrator);
		temperatureSetter.setIsothermalButtonsVisibility(false);
		temperatureSetter.setPrecision(4);
		temperatureSetter.setMinimum(0.0);
		temperatureSetter.setMaximum(1.0);
		temperatureSetter.setSliderMajorValues(5);
		temperatureSetter.setIsothermal();
		temperatureSetter.setTemperature(sim.temperature);
		
		temperatureSetter.setSliderPostAction(new IAction() {
            public void actionPerformed() {
                waveVectorSlider.setMaximum(numWV);
                waveVectorSlider.setIntegrator(sim.integrator);
                //Array
                if(sim.integrator.getWaveVectorNum() >= numWV){
                	waveVectorSlider.setMaximum(numWV);
                	sim.integrator.setWaveVectorNum(numWV-1);
                
                }
               
                /*
                 * Harmonic Free Energy
                 */
                double AHarmonic =0;
                int wvNum = (int)waveVectorSlider.getWaveVectorNum();
                int eValNum = (int)eValSlider.getEValNum();
                
                coeffs = sim.nm.getWaveVectorFactory().getCoefficients();
       
                if (sim.integrator.isOneWV()){
                	
                	if (sim.integrator.isOneEVal()){
                		AHarmonic=0;
                		AHarmonic = coeffs[wvNum] * Math.log(omega2[wvNum][eValNum]*coeffs[wvNum] /
                				(sim.integrator.temperature*Math.PI));
                		
                	} else {
                		AHarmonic=0;
                		for (int j=0; j<omega2[0].length; j++){
                			if (!Double.isInfinite(omega2[wvNum][j])) {
                				AHarmonic += coeffs[wvNum] * Math.log(omega2[wvNum][j]*coeffs[wvNum] /
                						(sim.integrator.temperature*Math.PI));
                        	}
                		}
                	}
                	
                }
                AHarm.E(AHarmonic);
                DataInfoDouble dataInfoA = new DataInfoDouble("AHarmonic", Energy.DIMENSION);
                displayAHarmonic.putDataInfo(dataInfoA);
                displayAHarmonic.putData(AHarm);
                displayAHarmonic.repaint();
                
                getController().getSimRestart().getDataResetAction().actionPerformed();
            	
		    }
		});
       
		
		// end of Temperature Slider
		
		/*
		 * N Cell Slider
		 */
		final DeviceCellNum3DSlider cellSlider = new DeviceCellNum3DSlider(sim.getController());
		cellSlider.setBox(sim.box);
		cellSlider.setSpecies(sim.species);
		cellSlider.setMaximum(5);
		
        cellSlider.setNSliderPostAction(new IAction() {
        	
       	  	public void actionPerformed() {
       	  		
       	  		int n = (int)cellSlider.getNCellNum(); 
       	  	 
                if (oldn != n ) {
                	int [] nCells = new int[]{n,n,n};
                	sim.nm.setNCellNum(n);
                	
                	Boundary boundary = new BoundaryRectangularPeriodic(sim.getSpace(), n*sim.L);
                	sim.box.setBoundary(boundary);
                	sim.setNCells(nCells);
                	sim.coordinateDefinition.initializeCoordinates(nCells);
                	sim.waveVectorFactory.makeWaveVectors(sim.box);
                	sim.truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
                	sim.pTruncated.setTruncationRadius(sim.truncationRadius);
                	sim.latticeEnergy = sim.meterPE.getDataAsScalar();
                	                	
                	numWV = sim.nm.getOmegaSquared().length;
                    numEval = sim.nm.getOmegaSquared()[0].length;
                    
                    omega2 = new double[numWV][numEval];
                    eigenVectors = new double[numWV][numEval][numEval];
                    
                }            
                
                oldn = n;
                
                int wvNumUsed = (int)waveVectorSlider.getWaveVectorNum();
                
                sim.integrator.reset();
                
             	if (sim.integrator.waveVectorNum >= numWV){
            		waveVectorSlider.setMaximum(numWV);
            		sim.integrator.setWaveVectorNum(numWV-1);
            		wvNumUsed = sim.integrator.getWaveVectorNum();
            		
            	}
             	
                sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(), sim.waveVectorFactory.getCoefficients());
                sim.integrator.setEigenVectors(sim.nm.getEigenvectors());
                sim.integrator.setWaveVectors(sim.waveVectorFactory.getWaveVectors());
                sim.integrator.setWaveVectorCoefficients(sim.waveVectorFactory.getCoefficients());
           
            	eValSlider.setMaximum(numEval);
            	waveVectorSlider.setMaximum(numWV);
            	
       
        		//change for wave vectors
                
                wavevectorx = new double[numWV];
                wavevectory = new double[numWV];
                wavevectorz = new double[numWV];
                stringiWV = new String[numWV];
                
                if (sim.integrator.isOneWV()){
                	sim.integrator.setWaveVectorNum(wvNumUsed);
                	           	
                	for (int nWV=0; nWV<numWV; nWV++){
               			wavevectorx[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(0);
               			wavevectory[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(1);
               			wavevectorz[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(2);
               			
               			if(nWV==wvNumUsed){
               				stringiWV[nWV] = String.valueOf("< WV "+nWV+">");
           
                		} else {
                			stringiWV[nWV] = String.valueOf("  WV "+nWV);
                		}
                		
                	}
                } else {
                	
                   	for (int nWV=0; nWV<numWV; nWV++){
                		wavevectorx[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(0);
               			wavevectory[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(1);
               			wavevectorz[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(2);
                		
               			stringiWV[nWV] = String.valueOf("< WV "+nWV+">");
               	  	}
                }
                    
                dataWV = new DataDoubleArray[3];
                dataWV[0] = new DataDoubleArray(new int[]{numWV}, wavevectorx);
                dataWV[1] = new DataDoubleArray(new int[]{numWV}, wavevectory);
                dataWV[2] = new DataDoubleArray(new int[]{numWV}, wavevectorz);
                wvTable = new DataTable(dataWV);
                
                DataInfoDoubleArray columnInfoWVx = new DataInfoDoubleArray("x", Null.DIMENSION, new int[]{numWV});
                DataInfoDoubleArray columnInfoWVy = new DataInfoDoubleArray("y", Null.DIMENSION, new int[]{numWV});
                DataInfoDoubleArray columnInfoWVz = new DataInfoDoubleArray("z", Null.DIMENSION, new int[]{numWV});
              
                DataInfo dataInfoTableWV = new DataInfoTable("Wave Vector", 
                		new DataInfoDoubleArray[]{columnInfoWVx,columnInfoWVy,columnInfoWVz}, numWV, stringiWV);
                //sinkWV.putDataInfo(dataInfoTableWV);
                sinkWV.putData(wvTable);
                
                // end of change for wave vectors
                
                waveVectorPostAction.actionPerformed();
                eValPostAction.actionPerformed();
                
                getController().getSimRestart().getDataResetAction().actionPerformed();
                getDisplayBox(sim.box).repaint();
       	  	}
       	  	
       	  	int oldn = (int)cellSlider.getNCellNum();
        });
        
        //End of N-Cell Slider
		
        /*
         *  Wave vectors Slider
         */ 

        waveVectorSlider = new DeviceWaveVectorSlider(sim.getController());
        waveVectorSlider.setMinimum(0);
        waveVectorSlider.setMaximum(sim.nm.getOmegaSquared().length);
        //waveVectorSlider.setOneWVButtonsVisibility(false);
        waveVectorSlider.setIntegrator(sim.integrator);
        
        waveVectorPostAction = new IAction() {
        	
        	public void actionPerformed() {
        		
                int wvNumUsed = (int)waveVectorSlider.getWaveVectorNum();
                
        		//change for wave vectors
                numWV = sim.nm.getOmegaSquared().length;
                numEval = sim.nm.getOmegaSquared()[0].length;
                
                wavevectorx = new double[numWV];
                wavevectory = new double[numWV];
                wavevectorz = new double[numWV];
                stringiWV = new String[numWV];
                
                if (sim.integrator.isOneWV()){
                	
                	for (int nWV=0; nWV<numWV; nWV++){
               			wavevectorx[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(0);
               			wavevectory[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(1);
               			wavevectorz[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(2);
               			
               			if(nWV==wvNumUsed){
               				stringiWV[nWV] = String.valueOf("< WV "+nWV+">");
           
                		} else {
                			stringiWV[nWV] = String.valueOf("  WV "+nWV);
                		}
                		
                	}
                } else {
                   	for (int nWV=0; nWV<numWV; nWV++){
                		wavevectorx[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(0);
               			wavevectory[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(1);
               			wavevectorz[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(2);
                		
               			stringiWV[nWV] = String.valueOf("< WV "+nWV+">");
               	  	}
                }
                        
                dataWV[0] = new DataDoubleArray(new int[]{numWV}, wavevectorx);
                dataWV[1] = new DataDoubleArray(new int[]{numWV}, wavevectory);
                dataWV[2] = new DataDoubleArray(new int[]{numWV}, wavevectorz);
                wvTable = new DataTable(dataWV);
                
                DataInfoDoubleArray columnInfoWVx = new DataInfoDoubleArray("x", Null.DIMENSION, new int[]{numWV});
                DataInfoDoubleArray columnInfoWVy = new DataInfoDoubleArray("y", Null.DIMENSION, new int[]{numWV});
                DataInfoDoubleArray columnInfoWVz = new DataInfoDoubleArray("z", Null.DIMENSION, new int[]{numWV});
                
                DataInfo dataInfoTableWV = new DataInfoTable("Wave Vector", 
                		new DataInfoDoubleArray[]{columnInfoWVx,columnInfoWVy,columnInfoWVz}, numWV, stringiWV);
                sinkWV.putDataInfo(dataInfoTableWV);
                sinkWV.putData(wvTable);
                // end of change for wave vectors
                
                eValPostAction.actionPerformed();
                
                getController().getSimRestart().getDataResetAction().actionPerformed();
                
		    }
		};
		
		ActionListener isOneWVListener = new ActionListener(){
			public void actionPerformed (ActionEvent event){
				waveVectorPostAction.actionPerformed();
				
			}
		};
        
		waveVectorSlider.setSliderPostAction(waveVectorPostAction);
		waveVectorSlider.addRadioGroupActionListener(isOneWVListener);

        // end wave vectors slider
        
        /*
         * Eigenvalues Slider
         */
        eValSlider = new DeviceEigenvaluesSlider(sim.getController());
        eValSlider.setMinimum(0);
        eValSlider.setMaximum(sim.nm.getOmegaSquared()[0].length);
        eValSlider.setIntegrator(sim.integrator);
        sim.integrator.setOneWV(true);
        
        eValPostAction = new IAction(){
        	public void actionPerformed() {
        		
        		int numWV = sim.nm.getOmegaSquared().length; 
        		int eValNum = (int)eValSlider.getEValNum();
                int wvNum = (int)waveVectorSlider.getWaveVectorNum();
                
                // eigenvector
                
                omega2 = new double[numWV][numEval];
                eVec = new double[numEval*numEval];
                so2 = new String[numEval*numEval];
                stringOmega2 = new String[numEval];
                
                for (int nEval1=0; nEval1<numEval; nEval1++){
                		omega2[wvNum][nEval1] = sim.nm.getOmegaSquared()[wvNum][nEval1];
                }
                
                
                if (eValSlider.isOneEVal()){
                	sim.integrator.reset();
                	sim.integrator.setOneEVal(true);
                    sim.integrator.setWaveVectorNum(wvNum);
                    sim.integrator.setEValNum(eValNum);
                    sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(), sim.waveVectorFactory.getCoefficients());
                    sim.integrator.setEigenVectors(sim.nm.getEigenvectors());
                	
                	
                	for (int nEval=0; nEval<numEval; nEval++){
                		for (int nEval2=0; nEval2<numEval; nEval2++){
                			if (nEval==eValNum){
                				stringOmega2[nEval] = "<"+String.valueOf(omega2[wvNum][nEval])+">";
               		
                			} else {
                				stringOmega2[nEval] = String.valueOf(omega2[wvNum][nEval]);
                			}
                			
                			if(nEval2==0){
                				so2[(nEval*numEval)+nEval2] = stringOmega2[nEval];
                			} else {
                				so2[(nEval*numEval)+nEval2] = " ";
                			}
                			
                			eVec[(nEval*numEval)+nEval2] = sim.nm.eigenvectors[sim.integrator.getWaveVectorNum()][nEval][nEval2];
                		}                   		
                	}
                	
                	meterHarmonicCoordinate.setEigenvectors(eigenVectors[wvNum][eValNum]);
                    meterHarmonicCoordinate.setWaveVector(sim.waveVectorFactory.getWaveVectors()[wvNum]);
                	
                } else {
                	sim.integrator.reset();
                	sim.integrator.setOneWV(true);
                	sim.integrator.setOneEVal(false);
                	sim.integrator.setWaveVectorNum(wvNum);
                    sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(), sim.waveVectorFactory.getCoefficients());
                    sim.integrator.setEigenVectors(sim.nm.getEigenvectors());
                    
                	for (int nEval=0; nEval<numEval; nEval++){
                		for (int nEval2=0; nEval2<numEval; nEval2++){
                			
                			stringOmega2[nEval] = "<"+String.valueOf(omega2[wvNum][nEval])+">";
                			                			
                			if(nEval2==0){
                				so2[(nEval*numEval)+nEval2] = stringOmega2[nEval];
                			} else {
                				so2[(nEval*numEval)+nEval2] = " ";
                			}
                			
                			eVec[(nEval*numEval)+nEval2] = sim.nm.eigenvectors[sim.integrator.getWaveVectorNum()][nEval][nEval2];
                		}
                		
                	}
                }
                
                
                dataEVec = new DataDoubleArray[1];
                dataEVec[0] = new DataDoubleArray(new int[]{numEval*numEval},eVec);
                eigenTable = new DataTable(dataEVec);
                          
                DataInfoDoubleArray columnInfoEigen = new DataInfoDoubleArray("Eigenvector", Null.DIMENSION, new int[]{numEval*numEval});
                DataInfo dataInfoTableEigen = new DataInfoTable("Eigenvector", new DataInfoDoubleArray[]{columnInfoEigen}, (numEval*numEval), so2);
                sinkEigen.putDataInfo(dataInfoTableEigen);
                sinkEigen.putData(eigenTable);
                
                // end of eigenvectors
                
                /*
                 * Harmonic Free Energy
                 */
                double AHarmonic =0;
                
                coeffs = sim.nm.getWaveVectorFactory().getCoefficients();
       
                if (sim.integrator.isOneWV()){
                	
                	if (sim.integrator.isOneEVal()){
                		AHarmonic=0;
                		AHarmonic = coeffs[wvNum] * Math.log(omega2[wvNum][eValNum]*coeffs[wvNum] /
                				(sim.integrator.temperature*Math.PI));
                		
                	} else {
                		AHarmonic=0;
                		for (int j=0; j<omega2[0].length; j++){
                			if (!Double.isInfinite(omega2[wvNum][j])) {
                				AHarmonic += coeffs[wvNum] * Math.log(omega2[wvNum][j]*coeffs[wvNum] /
                						(sim.integrator.temperature*Math.PI));
                        	}
                		}
                	}
                	
                }
                AHarm.E(AHarmonic);
                DataInfoDouble dataInfoA = new DataInfoDouble("AHarmonic", Energy.DIMENSION);
                displayAHarmonic.putDataInfo(dataInfoA);
                displayAHarmonic.putData(AHarm);
                displayAHarmonic.repaint();
                
                
                getController().getSimRestart().getDataResetAction().actionPerformed();
		    }
		};	
		
		ActionListener isOneEvalListener = new ActionListener(){
			public void actionPerformed (ActionEvent event){
				eValPostAction.actionPerformed();
			}
		};
        
		eValSlider.setSliderPostAction(eValPostAction);
		eValSlider.addRadioGroupActionListener(isOneEvalListener);
		
        
        //End of Eigenvalues Slider
        

        
		/*
		 * Display Table WV-Eigenvalues and Eigenvalues-Eigenvectors
		 */
		
		//tabbed-pane for wave vectors
        numWV = sim.nm.getOmegaSquared().length;
        numEval = sim.nm.getOmegaSquared()[0].length;
        
        wavevectorx = new double[numWV];
        wavevectory = new double[numWV];
        wavevectorz = new double[numWV];
        
        stringiWV = new String[numWV];
        
        if (waveVectorSlider.isOneWV()){
        	int wvNumUsed = sim.integrator.getWaveVectorNum();
        	
        	for (int nWV=0; nWV<numWV; nWV++){
       			wavevectorx[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(0);
       			wavevectory[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(1);
       			wavevectorz[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(2);
       			
       			if(nWV==wvNumUsed){
       				stringiWV[nWV] = String.valueOf("< WV "+nWV+">");
   
        		} else {
        			stringiWV[nWV] = String.valueOf("  WV "+nWV);
        		}
        		
        	}
        } else {
           	for (int nWV=0; nWV<numWV; nWV++){
        		wavevectorx[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(0);
       			wavevectory[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(1);
       			wavevectorz[nWV] = sim.waveVectorFactory.getWaveVectors()[nWV].getX(2);
        		
       			stringiWV[nWV] = String.valueOf("< WV "+nWV+">");
       	  	}
        }
 
        displayTableWV = new DisplayTable();
        displayTableWV.setTransposed(false);

        dataWV = new DataDoubleArray[3];
        dataWV[0] = new DataDoubleArray(new int[]{numWV}, wavevectorx);
        dataWV[1] = new DataDoubleArray(new int[]{numWV}, wavevectory);
        dataWV[2] = new DataDoubleArray(new int[]{numWV}, wavevectorz);
        wvTable = new DataTable(dataWV);
        
        DataInfoDoubleArray columnInfoWVx = new DataInfoDoubleArray("x", Null.DIMENSION, new int[]{numWV});
        DataInfoDoubleArray columnInfoWVy = new DataInfoDoubleArray("y", Null.DIMENSION, new int[]{numWV});
        DataInfoDoubleArray columnInfoWVz = new DataInfoDoubleArray("z", Null.DIMENSION, new int[]{numWV});
        
        DataInfo dataInfoTableWV = new DataInfoTable("Wave Vector", 
        		new DataInfoDoubleArray[]{columnInfoWVx,columnInfoWVy,columnInfoWVz}, numWV, stringiWV);
        sinkWV = displayTableWV.getDataTable().makeDataSink(dataInfoTableWV);
        sinkWV.putDataInfo(dataInfoTableWV);
        sinkWV.putData(wvTable);
        
        getPanel().tabbedPane.add("Wave Vector", displayTableWV.graphic());
        
        /*
         * tabbed-pane for eigenvectors with corresponding omega2
         */
        // eigenvector
        omega2 = new double[numWV][numEval];
        eigenVectors = new double[numWV][numEval][numEval];
        eVec = new double[numEval*numEval];
        so2 = new String[numEval*numEval];
        stringOmega2 = new String[numEval];
        
        for (int nWV=0; nWV<numWV; nWV++){
        	for (int nEval1=0; nEval1<numEval; nEval1++){
        		omega2[nWV][nEval1] = sim.nm.getOmegaSquared()[nWV][nEval1];
        		
        		for (int nEval2=0; nEval2<numEval; nEval2++){
        			eigenVectors[nWV][nEval1][nEval2] = sim.nm.eigenvectors[nWV][nEval1][nEval2];
        		}
        	}
        }
        // one eigenvalue check
        if (eValSlider.isOneEVal()){
        	int eValNumUsed = sim.integrator.getEValNum();
        	
        	for (int nEval=0; nEval<numEval; nEval++){
        		for (int nEval2=0; nEval2<numEval; nEval2++){
        			if (nEval==eValNumUsed){
        				stringOmega2[nEval] = "<"+String.valueOf(omega2[sim.integrator.getWaveVectorNum()][nEval])+">";
        		
        			} else {
        				stringOmega2[nEval] = String.valueOf(omega2[sim.integrator.getWaveVectorNum()][nEval]);
        			}
        			
        			if(nEval2==0){
        				so2[(nEval*numEval)+nEval2] = stringOmega2[nEval];
            		} else {
            			so2[(nEval*numEval)+nEval2] = " ";
            		}
        			
        			eVec[(nEval*numEval)+nEval2] = eigenVectors[sim.integrator.getWaveVectorNum()][nEval][nEval2];
        		}
        		
        	}
        	        	
        } else {
        	for (int nEval=0; nEval<numEval; nEval++){
        		for (int nEval2=0; nEval2<numEval; nEval2++){
        			stringOmega2[nEval] = "<"+String.valueOf(omega2[sim.integrator.getWaveVectorNum()][nEval])+">";
        		        			
        			if(nEval2==0){
        				so2[(nEval*numEval)+nEval2] = stringOmega2[nEval];
            		} else {
            			so2[(nEval*numEval)+nEval2] = " ";
            		}
        			
        			eVec[(nEval*numEval)+nEval2] = eigenVectors[sim.integrator.getWaveVectorNum()][nEval][nEval2];
        		}
        	}
        }
        
      
        displayTableEigen = new DisplayTable();
        displayTableEigen.setTransposed(false);
        
        dataEVec = new DataDoubleArray[1];
        dataEVec[0] = new DataDoubleArray(new int[]{numEval*numEval},eVec);
        eigenTable = new DataTable(dataEVec);
                  
        DataInfoDoubleArray columnInfoEigen = new DataInfoDoubleArray("Eigenvector", Null.DIMENSION, new int[]{numEval*numEval});
        DataInfo dataInfoTableEigen = new DataInfoTable("Eigenvector", new DataInfoDoubleArray[]{columnInfoEigen}, (numEval*numEval), so2);
        sinkEigen = displayTableEigen.getDataTable().makeDataSink(dataInfoTableEigen);
        sinkEigen.putDataInfo(dataInfoTableEigen);
        sinkEigen.putData(eigenTable);
        
        getPanel().tabbedPane.add("Omega^2", displayTableEigen.graphic());
        // end of eigenvectors
        
        
        
        /*
         * Harmonic Free Energy
         */
       
        AHarm = new DataDouble();
        double AHarmonic =0;
        
        coeffs = sim.nm.getWaveVectorFactory().getCoefficients();
        
        for(int i=0; i<omega2.length; i++) {
        	for (int j=0; j<omega2[0].length; j++){
                if (!Double.isInfinite(omega2[i][j])) {
                    AHarmonic += coeffs[i] * Math.log(omega2[i][j]*coeffs[i] /
                            (sim.integrator.temperature*Math.PI));
                }
        	}
        }
        
        if (sim.integrator.isOneWV()){
        	AHarmonic=0;
        	int wvNum = (int)waveVectorSlider.getWaveVectorNum();
        	for (int j=0; j<omega2[0].length; j++){
                if (!Double.isInfinite(omega2[wvNum][j])) {
                    AHarmonic += coeffs[wvNum] * Math.log(omega2[wvNum][j]*coeffs[wvNum] /
                            (sim.integrator.temperature*Math.PI));
                }
        	}
        	
        }
        AHarm.E(AHarmonic);
        
        displayAHarmonic = new DisplayTextBox();
        displayAHarmonic.setPrecision(10);
        dataInfoA = new DataInfoDouble("AHarmonic", Energy.DIMENSION);
        displayAHarmonic.putDataInfo(dataInfoA);
        displayAHarmonic.putData(AHarm);
        displayAHarmonic.setLabel("Harmonic Free Energy");
              
        getController().getDataStreamPumps().add(hePump);
        
        resetAction = new IAction(){
        	public void actionPerformed(){
        		heDisplay.putData(heAccumulator.getData());
        		heDisplay.repaint();
        		getDisplayBox(sim.box).graphic().repaint();
        	}
        };
        
        meterHarmonicCoordinate.setEigenvectors(eigenVectors[sim.integrator.getWaveVectorNum()][sim.integrator.eValNum]);
        meterHarmonicCoordinate.setWaveVector(sim.waveVectorFactory.getWaveVectors()[sim.integrator.getWaveVectorNum()]);
        
        this.getController().getReinitButton().setPostAction(resetAction);
        this.getController().getResetAveragesButton().setPostAction(resetAction);
        
        add(cellSlider);
        add(temperatureSetter);
        add(waveVectorSlider);
        add(displayAHarmonic);
        add(eValSlider);
        add(ePlot);
        add(eQPlot);
        add(heDisplay);
        
	}
	
	public static void main(String[] args){
		Space sp = Space.getInstance(3);
		NormalModeAnalysisDisplay3DGraphic simGraphic = new NormalModeAnalysisDisplay3DGraphic(new NormalModeAnalysisDisplay3D(sp), sp);
		SimulationGraphic.makeAndDisplayFrame(simGraphic.getPanel(), APP_NAME);
		
	}
	
	public static class Applet extends javax.swing.JApplet {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public void init(){
			getRootPane().putClientProperty(APP_NAME, Boolean.TRUE);
			Space sp = Space.getInstance(3);
			NormalModeAnalysisDisplay3DGraphic nm3Dgraphic = new NormalModeAnalysisDisplay3DGraphic(new NormalModeAnalysisDisplay3D(sp), sp);
			getContentPane().add(nm3Dgraphic.getPanel());
		}
	}
	
	private DeviceThermoSlider temperatureSetter; 
	private DeviceWaveVectorSlider waveVectorSlider;
	private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "3-D Harmonic Oscillator";
	private static final int REPAINT_INTERVAL = 10;
	protected NormalModeAnalysisDisplay3D sim;
	
	
	protected DataTable wvTable, eigenTable;
	protected DataDoubleArray[] dataWV, dataEVec;
	protected final DisplayTable displayTableWV, displayTableEigen;
	protected IDataSink sinkWV, sinkEigen;
	protected String[] stringiWV, stringOmega2, so2;
	protected double[] coeffs;
	protected IData AHarm;
	protected DisplayTextBox displayAHarmonic;
	protected DataInfoDouble dataInfoA; 
	protected DeviceEigenvaluesSlider eValSlider;
	protected int numWV, numEval;
	protected double[] wavevectorx, wavevectory, wavevectorz, eVec;
	
	protected double[][] omega2;
	protected double[][][] eigenVectors;
	protected final IAction waveVectorPostAction, eValPostAction, resetAction;
	protected MeterHarmonicCoordinate meterHarmonicCoordinate;
	
}
