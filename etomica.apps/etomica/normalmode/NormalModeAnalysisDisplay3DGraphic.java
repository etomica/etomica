package etomica.normalmode;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import etomica.api.IAction;
import etomica.api.IData;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataInfo;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.data.DataTag;
import etomica.data.IDataSink;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTable;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataTable.DataInfoTable;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DeviceThermoSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayTextBoxesCAE;
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
public class NormalModeAnalysisDisplay3DGraphic extends SimulationGraphic {

	public NormalModeAnalysisDisplay3DGraphic(final NormalModeAnalysisDisplay3D simulation, Space space) {
		
		super(simulation, TABBED_PANE, APP_NAME,REPAINT_INTERVAL, space, simulation.getController());
		this.sim = simulation;
		
		DataSourceCountTime timeCounter = new DataSourceCountTime(sim.integrator);
		
		
		/*
		 * harmonic energy                                                                            
		 */
		MeterHarmonicEnergy heMeter = new MeterHarmonicEnergy(sim.coordinateDefinition, sim.nm);
		heMeter.setBox(sim.box);
		
		AccumulatorHistory heHistory = new AccumulatorHistory();
		heHistory.setTimeDataSource(timeCounter);
	    
		final AccumulatorAverageCollapsing heAccumulator = new AccumulatorAverageCollapsing();
        heAccumulator.setPushInterval(10);
        DataFork heFork = new DataFork(new IDataSink[]{heHistory, heAccumulator});
        DataPump hePump = new DataPump(heMeter, heFork);
        sim.integrator.addIntervalAction(hePump);
        sim.integrator.setActionInterval(hePump, 60);
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
		temperatureSetter = new DeviceThermoSlider(sim.getController());
		temperatureSetter.setIsothermalButtonsVisibility(false);
		temperatureSetter.setPrecision(2);
		temperatureSetter.setMinimum(0.0);
		temperatureSetter.setMaximum(1.0);
		temperatureSetter.setSliderMajorValues(5);
		temperatureSetter.setIntegrator(sim.integrator);
		temperatureSetter.setIsothermal();
		temperatureSetter.setTemperature(sim.temperature);
		
		temperatureSetter.setSliderPostAction(new IAction() {
            public void actionPerformed() {
            	                
                
                int numWV = sim.nm.getOmegaSquared(null).length;
                int numEval = sim.nm.getOmegaSquared(null)[0].length;
                int totalDim = numWV*numEval;
        		
                double[][] omega2 = new double[numWV][numEval];
        		double[] o2 = new double[totalDim];
        		String[] sWV = new String[totalDim]; 
                stringWV = new String[numWV][numEval];

                waveVectorSlider.setMaximum(numWV);
                waveVectorSlider.setIntegrator(sim.integrator);
                //Array
                if(sim.integrator.getWaveVectorNum() >= numWV){
                	waveVectorSlider.setMaximum(numWV);
                	sim.integrator.setWaveVectorNum(numWV-1);
                }
                /*
                if (sim.integrator.isOneWV()){
                	int wvNumUsed = sim.integrator.getWaveVectorNum();
                	
                	for (int nWV=0; nWV<numWV; nWV++){
                		for (int nEval=0; nEval<numEval; nEval++){
                			omega2[nWV][nEval] = sim.nm.getOmegaSquared(sim.box)[nWV][nEval];

                			if(nWV==wvNumUsed){
                					
                				if(nEval==0){  // assign the WV value to stringWV[nWV][0]
                					stringWV[nWV][nEval]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(0))+", "
                					+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(1))+", "
                					+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(2))+">";
                				
                				} else {
                					stringWV[nWV][nEval] = " ";
                				}

                			} else {
                				if(nEval==0){
                					stringWV[nWV][nEval]=String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(0))+", "
                					+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(1))+", "
                					+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(2));
                				
                				} else {
                					stringWV[nWV][nEval] = " ";
                				
                				}
                			}
                			
                			o2[(nWV*numEval) + nEval] = omega2[nWV][nEval];
                			sWV[(nWV*numEval) + nEval] = stringWV[nWV][nEval];
                		}
                	}
                               
                } else {
                
                	for (int nWV=0; nWV<numWV; nWV++){
                		for (int nEval=0; nEval<numEval; nEval++){
                			
                			omega2[nWV][nEval] = sim.nm.getOmegaSquared(sim.box)[nWV][nEval];
                			if(nEval==0){
                				stringWV[nWV][nEval]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(0))+", "
                				+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(1))+", "
                				+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(2))+">";
                			
                			} else {
                				stringWV[nWV][nEval] = " ";
                			}
                			o2[(nWV*numEval) + nEval] = omega2[nWV][nEval];
                			sWV[(nWV*numEval) + nEval] = stringWV[nWV][nEval];
                		}
                	}
                }
         
                
                // tabbed-pane for omega2's with corresponding wave vectos
                 
                dataOmega2 = new DataDoubleArray[1];
                dataOmega2[0] = new DataDoubleArray(new int[]{totalDim},o2);
                wvTable = new DataTable(dataOmega2);
                 
                DataInfoDoubleArray columnInfoWV = new DataInfoDoubleArray("Omega^2", Null.DIMENSION, new int[]{totalDim});
                DataInfo dataInfoTableWV = new DataInfoTable("Omega^2", new DataInfoDoubleArray[]{columnInfoWV}, totalDim, sWV);
                sinkWV.putDataInfo(dataInfoTableWV);
                sinkWV.putData(wvTable);
                
                getPanel().tabbedPane.add("Wave Vector", displayTableWV.graphic());
                
                
                //  tabbed-pane for eigenvectors with corresponding omega2
                 
                
                
                double[][][] eigenVectors = new double[numWV][numEval][numEval];
                double[] eVec = new double[numEval*numEval];
                String[] so2 = new String[numEval*numEval];
                stringOmega2 = new String[numEval];
                
                for (int nWV=0; nWV<numWV; nWV++){
                	for (int nEval1=0; nEval1<numEval; nEval1++){
                		for (int nEval2=0; nEval2<numEval; nEval2++){
                			eigenVectors[nWV][nEval1][nEval2] = sim.nm.eigenvectors[nWV][nEval1][nEval2];
                		}
                	}
                }
                // one eigenvalue check
                if (sim.integrator.isOneEVal()){
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
                
                
                dataEVec = new DataDoubleArray[1];
                dataEVec[0] = new DataDoubleArray(new int[]{numEval*numEval},eVec);
                eigenTable = new DataTable(dataEVec);
        
                DataInfoDoubleArray columnInfoEigen = new DataInfoDoubleArray("Eigenvector", Null.DIMENSION, new int[]{numEval*numEval});
                DataInfo dataInfoTableEigen = new DataInfoTable("Eigenvector", new DataInfoDoubleArray[]{columnInfoEigen}, (numEval*numEval), so2);
                sinkEigen.putDataInfo(dataInfoTableEigen);
                sinkEigen.putData(eigenTable);
            	*/
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
                AHarm.E(AHarmonic);
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
		cellSlider.setMaximum(sim.getN());
		cellSlider.setNCellNum(sim.n);
		
        cellSlider.setNSliderPostAction(new IAction() {
        	
       	  	public void actionPerformed() {
       	  		
       	  		int n = (int)cellSlider.getNCellNum(); 
       	  	 
                if (oldn != n ) {
                	int [] nCells = new int[]{n,n,n};
           
                	Boundary boundary = new BoundaryRectangularPeriodic(sim.getSpace(), n*sim.L);
                	sim.box.setBoundary(boundary);
                	
                	sim.waveVectorFactory.makeWaveVectors(sim.box);
                	sim.coordinateDefinition.initializeCoordinates(nCells);
                	
                }            
                
                oldn = n;
                try {
                	sim.integrator.reset();
            	
                	sim.integrator.setWaveVectors(sim.waveVectorFactory.getWaveVectors());
                	sim.integrator.setWaveVectorCoefficients(sim.waveVectorFactory.getCoefficients());
                	sim.nm.setNCellNum(n);
                	sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(sim.box), sim.waveVectorFactory.getCoefficients());
                	sim.integrator.setEigenVectors(sim.nm.getEigenvectors(sim.box));
                	
                	eValSlider.setMaximum(sim.nm.getOmegaSquared(null)[0].length);
                    
                	if (eValSlider.isOneEVal()){
                		sim.integrator.setOneEVal(eValSlider.isOneEVal());
                    	sim.integrator.setEValNum((int)eValSlider.getEValNum());
                       	sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(sim.box), sim.waveVectorFactory.getCoefficients());
                    	sim.integrator.setEigenVectors(sim.nm.getEigenvectors(sim.box));
                    }
                	
                	int m = sim.nm.getOmegaSquared(null).length;
                	double[] omega2 = new double[m];
                
                	waveVectorSlider.setMaximum(m);
                	waveVectorSlider.setIntegrator(sim.integrator);
                	//Array
                	if(sim.integrator.getWaveVectorNum() >= m){
                		waveVectorSlider.setMaximum(m);
                		sim.integrator.setWaveVectorNum(m-1);
                	}
                /*
                   	if (sim.integrator.isOneWV()){
                		int wvNumUsed = sim.integrator.getWaveVectorNum();
                		for (int i=0; i<m; i++){
                			omega2[i] = sim.nm.getOmegaSquared(sim.box)[i][0];
                    	
                			if(i==wvNumUsed){
                				stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(0))+", "
                				+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(1))+", "
                				+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(2))+">";
                			} else {
                				stringWV[i]= String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(0))+", "
                				+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(1))+", "
                				+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(2));
                			}
                		}
                	} else {
                
                		for (int i=0; i<m; i++){
                			omega2[i] = sim.nm.getOmegaSquared(sim.box)[i][0];
                			stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(0))+", "
                			+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(1))+", "
                			+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(2))+">";
                		}
                	}
                	data[0] = new DataDoubleArray(new int[]{m},omega2);
                	omega2Table = new DataTable(data);
               
                  	DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Omega^2", Null.DIMENSION, new int[]{m});
                	DataInfo dataInfo = new DataInfoTable("Omega^2", new DataInfoDoubleArray[]{columnInfo}, m, stringWV);
                	sink.putDataInfo(dataInfo);
                	sink.putData(omega2Table);
                	
                	
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
                
                }
                
                catch (ConfigurationOverlapException e) {
                    throw new RuntimeException(e);
                }
                
                getController().getSimRestart().getDataResetAction().actionPerformed();
                getDisplayBox(sim.box).repaint();
       	  	}
       	  	
       	  	int oldn = (int)cellSlider.getNCellNum();
        });
        
        //End of N-Cell Slider
		
        /*
         * Eigenvalues Slider
         */
        eValSlider = new DeviceEigenvaluesSlider(sim.getController());
        eValSlider.setMinimum(0);
        eValSlider.setMaximum(sim.nm.getOmegaSquared(null)[0].length);
        eValSlider.setIntegrator(sim.integrator);
        sim.integrator.setOneWV(true);
        final IAction eValPostAction = new IAction(){
        	public void actionPerformed() {
        		
                int eValNum = (int)eValSlider.getEValNum();
                if (eValSlider.isOneEVal()){
                	sim.integrator.setOneEVal(eValSlider.isOneEVal());
                	
                	if (eValNum != eValOld){
                		sim.integrator.reset();
                		sim.integrator.setEValNum(eValNum);
                       	sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(sim.box), sim.waveVectorFactory.getCoefficients());
                    	sim.integrator.setEigenVectors(sim.nm.getEigenvectors(sim.box));
                	}
                	
                	eValOld = eValNum;
                	
                	sim.integrator.reset();
            		sim.integrator.setEValNum(eValNum);
                   	sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(sim.box), sim.waveVectorFactory.getCoefficients());
                	sim.integrator.setEigenVectors(sim.nm.getEigenvectors(sim.box));
                	
	            } else {
	            	sim.integrator.reset();
	            	sim.integrator.setOneEVal(false);
	              	sim.integrator.setOmegaSquared(sim.nm.getOmegaSquared(sim.box), sim.waveVectorFactory.getCoefficients());
                	sim.integrator.setEigenVectors(sim.nm.getEigenvectors(sim.box));
	            }
                
            	                
                getController().getSimRestart().getDataResetAction().actionPerformed();
                
		    }
        	int eValOld = (int)eValSlider.getEValNum();
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
         *  Wave vectors Slider
         */ 

        waveVectorSlider = new DeviceWaveVectorSlider(sim.getController());
        waveVectorSlider.setMinimum(0);
        waveVectorSlider.setMaximum(sim.nm.getOmegaSquared(null).length-1);
        waveVectorSlider.setOneWVButtonsVisibility(false);
        waveVectorSlider.setIntegrator(sim.integrator);
        
        final IAction waveVectorPostAction = new IAction() {
        	
        	public void actionPerformed() {
        		
                int m = sim.nm.getOmegaSquared(sim.box).length;
                double[] omega2 = new double[m];
                
                int wvNumUsed = (int)waveVectorSlider.getWaveVectorNum();
                /*
                if (waveVectorSlider.isOneWV()){
	                for (int i=0; i<m; i++){
	                	
	                	omega2[i] = sim.nm.getOmegaSquared(sim.box)[i][0];
	                	if (i==wvNumUsed){
	                		stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(0))+", "
	                		+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(1))+", "
	                		+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(2))+">";
	                		
	                	} else {
	                	
	                		stringWV[i]=String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(0))+", "
	                		+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(1))+", "
	                		+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(2));
	                	}
	                }
                } else {
                	for (int i=0; i<m; i++){
                    	omega2[i] = sim.nm.getOmegaSquared(sim.box)[i][0];
                    	stringWV[i]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(0))+", "
                    	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(1))+", "
                    	+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[i].x(2))+">";	
                    }
                }
                
                
                data[0] = new DataDoubleArray(new int[]{m},omega2);
                omega2Table = new DataTable(data);
                
                DataInfoDoubleArray columnInfo = new DataInfoDoubleArray("Omega^2", Null.DIMENSION, new int[]{m});
                DataInfo dataInfo = new DataInfoTable("Omega^2", new DataInfoDoubleArray[]{columnInfo}, m, stringWV);
                sink.putDataInfo(dataInfo);
                sink.putData(omega2Table);
            	             */   
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
		 * Display Table WV-Eigenvalues and Eigenvalues-Eigenvectors
		 */
        int numWV = sim.nm.getOmegaSquared(null).length;
        int numEval = sim.nm.getOmegaSquared(null)[0].length;
        int totalDim = numWV*numEval;
		
        double[][] omega2 = new double[numWV][numEval];
		double[] o2 = new double[totalDim];
		String[] sWV = new String[totalDim]; 
        stringWV = new String[numWV][numEval];

        if (sim.integrator.isOneWV()){
        	int wvNumUsed = sim.integrator.getWaveVectorNum();
        	
        	for (int nWV=0; nWV<numWV; nWV++){
        		for (int nEval=0; nEval<numEval; nEval++){
        			omega2[nWV][nEval] = sim.nm.getOmegaSquared(sim.box)[nWV][nEval];

        			if(nWV==wvNumUsed){
        					
        				if(nEval==0){  // assign the WV value to stringWV[nWV][0]
        					stringWV[nWV][nEval]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(0))+", "
        					+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(1))+", "
        					+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(2))+">";
        				
        				} else {
        					stringWV[nWV][nEval] = " ";
        				}

        			} else {
        				if(nEval==0){
        					stringWV[nWV][nEval]=String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(0))+", "
        					+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(1))+", "
        					+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(2));
        				
        				} else {
        					stringWV[nWV][nEval] = " ";
        				
        				}
        			}
        			
        			o2[(nWV*numEval) + nEval] = omega2[nWV][nEval];
        			sWV[(nWV*numEval) + nEval] = stringWV[nWV][nEval];
        		}
        	}
                       
        } else {
        
        	for (int nWV=0; nWV<numWV; nWV++){
        		for (int nEval=0; nEval<numEval; nEval++){
        			
        			omega2[nWV][nEval] = sim.nm.getOmegaSquared(sim.box)[nWV][nEval];
        			if(nEval==0){
        				stringWV[nWV][nEval]="<"+String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(0))+", "
        				+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(1))+", "
        				+ String.valueOf(sim.waveVectorFactory.getWaveVectors()[nWV].x(2))+">";
        			
        			} else {
        				stringWV[nWV][nEval] = " ";
        			}
        			o2[(nWV*numEval) + nEval] = omega2[nWV][nEval];
        			sWV[(nWV*numEval) + nEval] = stringWV[nWV][nEval];
        		}
        	}
        }
 
        /*
         * tabbed-pane for omega2's with corresponding wave vectos
         */
        dataOmega2 = new DataDoubleArray[1];
        dataOmega2[0] = new DataDoubleArray(new int[]{totalDim},o2);
        wvTable = new DataTable(dataOmega2);
         
        displayTableWV = new DisplayTable();
        sinkWV = displayTableWV.getDataTable().makeDataSink();
        
        displayTableWV.setTransposed(false);
        
        DataInfoDoubleArray columnInfoWV = new DataInfoDoubleArray("Omega^2", Null.DIMENSION, new int[]{totalDim});
        DataInfo dataInfoTableWV = new DataInfoTable("Omega^2", new DataInfoDoubleArray[]{columnInfoWV}, totalDim, sWV);
        sinkWV.putDataInfo(dataInfoTableWV);
        sinkWV.putData(wvTable);
        
        getPanel().tabbedPane.add("Wave Vector", displayTableWV.graphic());
        
        /*
         * tabbed-pane for eigenvectors with corresponding omega2
         */
        
        
        double[][][] eigenVectors = new double[numWV][numEval][numEval];
        double[] eVec = new double[numEval*numEval];
        String[] so2 = new String[numEval*numEval];
        stringOmega2 = new String[numEval];
        
        for (int nWV=0; nWV<numWV; nWV++){
        	for (int nEval1=0; nEval1<numEval; nEval1++){
        		for (int nEval2=0; nEval2<numEval; nEval2++){
        			eigenVectors[nWV][nEval1][nEval2] = sim.nm.eigenvectors[nWV][nEval1][nEval2];
        		}
        	}
        }
        // one eigenvalue check
        if (sim.integrator.isOneEVal()){
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
        //
        
        dataEVec = new DataDoubleArray[1];
        dataEVec[0] = new DataDoubleArray(new int[]{numEval*numEval},eVec);
        eigenTable = new DataTable(dataEVec);
        
        displayTableEigen = new DisplayTable();
        sinkEigen = displayTableEigen.getDataTable().makeDataSink();
        
        displayTableEigen.setTransposed(false);
        
        DataInfoDoubleArray columnInfoEigen = new DataInfoDoubleArray("Eigenvector", Null.DIMENSION, new int[]{numEval*numEval});
        DataInfo dataInfoTableEigen = new DataInfoTable("Eigenvector", new DataInfoDoubleArray[]{columnInfoEigen}, (numEval*numEval), so2);
        sinkEigen.putDataInfo(dataInfoTableEigen);
        sinkEigen.putData(eigenTable);
        
        getPanel().tabbedPane.add("Omega^2", displayTableEigen.graphic());
        
        //
        
        
        
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
        AHarm.E(AHarmonic);
        
        
        displayAHarmonic = new DisplayTextBox();
        displayAHarmonic.setPrecision(10);
        dataInfoA = new DataInfoDouble("AHarmonic", Energy.DIMENSION);
        displayAHarmonic.putDataInfo(dataInfoA);
        displayAHarmonic.putData(AHarm);
        displayAHarmonic.setLabel("Harmonic Free Energy");
              
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
        add(eValSlider);
        add(ePlot);
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
			NormalModeAnalysisDisplay3DGraphic nm1Dgraphic = new NormalModeAnalysisDisplay3DGraphic(new NormalModeAnalysisDisplay3D(sp), sp);
			getContentPane().add(nm1Dgraphic.getPanel());
		}
	}
	
	private DeviceThermoSlider temperatureSetter; 
	private DeviceWaveVectorSlider waveVectorSlider;
	private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "3-D Harmonic Oscillator";
	private static final int REPAINT_INTERVAL = 10;
	protected NormalModeAnalysisDisplay3D sim;
	
	
	protected DataTable wvTable, eigenTable;
	protected DataDoubleArray[] dataOmega2, dataEVec;
	protected final DisplayTable displayTableWV, displayTableEigen;
	protected IDataSink sinkWV, sinkEigen;
	protected String[][] stringWV;
	protected String[] stringOmega2;
	protected double[] coeffs;
	protected IData AHarm;
	protected DisplayTextBox displayAHarmonic;
	protected DataInfoDouble dataInfoA; 
	protected DeviceEigenvaluesSlider eValSlider;

}
