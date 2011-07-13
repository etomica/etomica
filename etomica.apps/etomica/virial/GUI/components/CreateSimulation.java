package etomica.virial.GUI.components;

import java.awt.Color;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiIntergroup;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverage;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P22CLJQ;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Spherical;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.CompoundDimension;
import etomica.units.CompoundUnit;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
import etomica.units.Kelvin;
import etomica.units.Liter;
import etomica.units.Mole;
import etomica.units.Pixel;
import etomica.units.Quantity;
import etomica.units.Unit;
import etomica.units.UnitRatio;
import etomica.units.Volume;
import etomica.util.Constants.CompassDirection;
import etomica.virial.ClusterAbstract;
import etomica.virial.MayerEGeneral;
import etomica.virial.MayerEHardSphere;
import etomica.virial.MayerESpherical;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesFactorySpheres;
import etomica.virial.SpeciesFactoryTangentSpheres;
import etomica.virial.GUI.containers.SimulationParameters;
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap;
import etomica.virial.simulations.SimulationVirialOverlapRejected;

public class CreateSimulation {
	
	
	private String CustomName = null;
	private SimulationEnvironment SimEnv;
	private Space space;
	
	private int[] MoleculeCount = new int[2];
	
	private SimulationGraphic simGraphic;
	private ParameterMapping[] potential = new ParameterMapping[2];
	
	private int nPoints;
	private boolean InterNonBondedPotentialFlag;
	private boolean Species1IntraNonBondedPotentialFlag;
	private boolean Species2IntraNonBondedPotentialFlag;
	private boolean Species1IntraBondedPotentialFlag;
	private boolean Species2IntraBondedPotentialFlag;
	
	private double temperature;
	private double SigmaHSRef;
	private int steps;
	private PotentialGroup pInterTargetGroup;
	private PotentialGroup pIntra1TargetGroup;
	private PotentialGroup pIntra2TargetGroup;
	
	Object[] Species1Potentials;
	Object[] Species2Potentials;
	private Constructor TempConstructor;
	
	private MayerGeneral fTarget;
	private MayerEGeneral eTarget;
	private MayerGeneralSpherical fTarget1; 
    private MayerESpherical eTarget1;
	
	
	private JFrame frame;
	
	public CreateSimulation(ParameterMapping[] Potential, int Molecule1Count, int Molecule2Count){
		
		this.potential[0] = Potential[0];
		if(Potential.length>1){
			this.potential[1] = Potential[1];
		}
		this.MoleculeCount[0]=Molecule1Count;
		this.MoleculeCount[1]=Molecule2Count;
		nPoints = Molecule1Count+Molecule2Count;
	}

	
	@SuppressWarnings("unchecked")
	public void runSimulation(SimulationEnvironment simenv) throws NoSuchMethodException{
		SimEnv = simenv;
		System.out.println(potential[0].getClass().getName());
		
		if(potential[1] != null){
			InterNonBondedPotentialFlag = true;
			if(potential[1].getPotentialSites().length > 1){
				if (etomica.api.IPotentialMolecular.class.isAssignableFrom(potential[1].getPotential())){
					Species2IntraNonBondedPotentialFlag = false;
					for(int i =0;i<potential[1].getPotentialSites().length;i++){
						if(potential[1].getPotentialSiteAtIndex(i).contains("CH3")){
							Species2IntraBondedPotentialFlag = true;
						}
						else{
							Species2IntraBondedPotentialFlag = false;
						}
					}
				}
				else{
					Species2IntraNonBondedPotentialFlag = true;
					for(int i =0;i<potential[1].getPotentialSites().length;i++){
						if(potential[1].getPotentialSiteAtIndex(i).contains("CH3")){
							Species2IntraBondedPotentialFlag = true;
						}
						else{
							Species2IntraBondedPotentialFlag = false;
						}
					}
				}
			}
		}
		else{
			InterNonBondedPotentialFlag = false;
			Species2IntraNonBondedPotentialFlag = false;
			Species2IntraBondedPotentialFlag = false;
		}
		
		
		if(potential[0].getPotentialSites().length > 1){
			if (etomica.api.IPotentialMolecular.class.isAssignableFrom(potential[0].getPotential())){
				Species1IntraNonBondedPotentialFlag = false;
				for(int i =0;i<potential[0].getPotentialSites().length;i++){
					if(potential[0].getPotentialSiteAtIndex(i).contains("CH3")){
						Species1IntraBondedPotentialFlag = true;
					}
					else{
						Species1IntraBondedPotentialFlag = false;
					}
				}
			}
			else{
				Species1IntraNonBondedPotentialFlag = true;
				for(int i =0;i<potential[0].getPotentialSites().length;i++){
					if(potential[0].getPotentialSiteAtIndex(i).contains("CH3")){
						Species1IntraBondedPotentialFlag = true;
					}
					else{
						Species1IntraBondedPotentialFlag = false;
					}
				}
			}
		}
		else{
			Species1IntraNonBondedPotentialFlag = false;
			Species1IntraBondedPotentialFlag = false;
		}
		
		
		
		//All Environment variables set first
		temperature = SimEnv.getTemperature();
		steps = SimEnv.getNoOfSteps();
		
		
		//SigmaHSRef will vary according to mixing rules
		SigmaHSRef = SimEnv.getSigmaHSRef();
		if(!InterNonBondedPotentialFlag){
			if(potential[0].getClass().getName().contains("P2Alkane")){
				if(Species1IntraBondedPotentialFlag){
					for(int i =0;i<potential[0].getPotentialSites().length;i++){
						if(potential[0].getPotentialSiteAtIndex(i).contains("CH2")){
							if(potential[0].getMoleculeDisplayName().contains("n-Alkane")){
								SigmaHSRef=potential[0].getDoubleDefaultParameters("SIGMACH3")+(SimEnv.getAlkane1Spheres()*0.5);
							}
							else{
								SigmaHSRef=potential[0].getDoubleDefaultParameters("SIGMACH3")+(0.5*3);}
						}
						else{
							
							SigmaHSRef=potential[0].getDoubleDefaultParameters("SIGMACH3")+(0.5*2);
						}
					}
				}
				else{
					if(potential[0].getPotentialSiteAtIndex(0).contains("CH3")){
						SigmaHSRef=potential[0].getDoubleDefaultParameters("SIGMACH3")+(0.5);
					}
					else{
						SigmaHSRef=potential[0].getDoubleDefaultParameters("SIGMACH4")+(0.5);
					}
				}
			}
		}
		//For Alkane mixtures
		
		 MayerHardSphere fRef = new MayerHardSphere(SigmaHSRef);
	     MayerEHardSphere eRef = new MayerEHardSphere(SigmaHSRef);
	     
	     //Set up the potentials
	     try{
	    	 if(!InterNonBondedPotentialFlag){
	    		 //Potential 1 for Intra Non-Bonded Potentials
	    		 if(Species1IntraNonBondedPotentialFlag){
	    			 pIntra1TargetGroup = new PotentialGroup(2);
	    			 
	    			 int NoOfIntraSites = potential[0].getPotentialSites().length;
	    			 int NoOfIntraSitePairs = Factorial(NoOfIntraSites)/(Factorial(NoOfIntraSites-2)*Factorial(2));
	    			 
	    			 Species1Potentials = new Object[potential[0].getPotentialSites().length + NoOfIntraSitePairs];
	    			
	    			 

	    			 @SuppressWarnings("rawtypes")
					 
	    			 
	    			 int numberofParameters = potential[0].getParametersArray().length;
	    			 Object[] ParamValueObj = new Object[numberofParameters+1];
	    			 
	    			 @SuppressWarnings("rawtypes")
	    			 Class[] ParamClass = new Class[numberofParameters+1];
	    			 
	    			 ParamClass[0] = ISpace.class;
					 for(int i = 1;i<=numberofParameters;i++){
						 ParamClass[i] = Double.TYPE;
					 }
	    			 
	    			 ParamValueObj[0] = potential[0].getSpace();
	    			 for(int index=0;index<NoOfIntraSites;index++){
	    				for(int j=0;j<numberofParameters;j++){
	    					ParamValueObj[j+1]=potential[0].getDoubleDefaultParameters(potential[0].getParametersArray()[j].toUpperCase()+potential[0].getPotentialSiteAtIndex(index));
						}
	    				if(potential[0].getPotential().getConstructors().length > 1){
	    					TempConstructor = potential[0].getPotential().getConstructor(ParamClass);
							Species1Potentials[index] = TempConstructor.newInstance(ParamValueObj);
	    				}
	    				else{
	    					TempConstructor = potential[0].getPotential().getConstructor(ISpace.class);
	    					Species1Potentials[index] = TempConstructor.newInstance(ParamValueObj[0]);
	    				}
	    			 }
	    			 
	    			 int[][] IntraSitePairs = new int[NoOfIntraSitePairs][2];
	    			 
	    			 int index = 0;
	    			 for(int l=0;l<NoOfIntraSites;l++){
	    				for(int k=NoOfIntraSites-1;k>l;k--){
	    					IntraSitePairs[index][0] = l;
	    					IntraSitePairs[index][1] = k;
	    					index++;
	    				}
	    			 }
	    			 for(int m=2;m<Species1Potentials.length;m++){
	    				 int speciesAIndex = IntraSitePairs[m-2][0];
	    				 int speciesBIndex = IntraSitePairs[m-2][1];
	    				 for(int j=0;j<numberofParameters;j++){
	    					String Param = potential[0].getParametersArray()[j];
	    					Double ParamValue = 0.0;
	    					Double ValueA;
	    					Double ValueB;
	    					ValueA = potential[0].getDoubleDefaultParameters(potential[0].getParametersArray()[j].toUpperCase()+potential[0].getPotentialSiteAtIndex(speciesAIndex));
	    					ValueB = potential[0].getDoubleDefaultParameters(potential[0].getParametersArray()[j].toUpperCase()+potential[0].getPotentialSiteAtIndex(speciesBIndex));
	    					if(Param.contains("SIGMA")){
	    						ParamValue = 0.5*(ValueA+ValueB);
	    					}
	    					if(Param.contains("EPSILON")){
	    						ParamValue = Math.sqrt(ValueA*ValueB);
	    					}
		    				ParamValueObj[j+1]=ParamValue;				
						 }
	    				 if(potential[0].getPotential().getConstructors().length > 1){
	    					 if(TempConstructor != null){
	    						 Species1Potentials[m] = TempConstructor.newInstance(ParamValueObj);
	    					 }
	    				 }
	    				 else{
	    					 if(TempConstructor != null){
	    						 
	    						 Species1Potentials[m] = TempConstructor.newInstance(ParamValueObj[0]);
	    					 }
	    				 }
	    			 }
	    			 fTarget = new MayerGeneral(pIntra1TargetGroup);
	    		     eTarget = new MayerEGeneral(pIntra1TargetGroup);

	    		 }
	    		 else{
	    			 	int numberofParameters = potential[0].getParametersArray().length;
						Object[] ParamValueObj = new Object[numberofParameters+1];
						
						@SuppressWarnings("rawtypes")
						Class[] ParamClass = new Class[numberofParameters+1];
						
						ParamValueObj[0] = potential[0].getSpace();
						for(int j=0;j<potential[0].getParametersArray().length;j++){
							ParamValueObj[j+1]=potential[0].getDoubleDefaultParameters(potential[0].getParametersArray()[j].toUpperCase()+potential[0].getPotentialSiteAtIndex(0));
						}
						
						
						ParamClass[0] = ISpace.class;
						for(int i = 1;i<=numberofParameters;i++){
							ParamClass[i] = Double.TYPE;
						}
						
	    			 	//All potentials whose constructors are 2 in number( one with simplece iSpace.class, and remaining with parameters!!)
	    			 	Species1Potentials = new Object[1];
						if(potential[0].getPotential().getConstructors().length > 1){
							Species1Potentials[0] = potential[0].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
						}
						else{
							//All potentials whose constructors simply consist of passing ISpaceClass
							Species1Potentials[0] = potential[0].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
						}
						if (etomica.api.IPotentialMolecular.class.isAssignableFrom(potential[0].getPotential())){
							fTarget = new MayerGeneral((IPotentialMolecular) Species1Potentials[0]);
			    		    eTarget = new MayerEGeneral((IPotentialMolecular)Species1Potentials[0]);
						}
						else{
							//pIntra1TargetGroup = new PotentialGroup(2);
							//pIntra1TargetGroup.addPotential((IPotentialAtomic)Species1Potentials[0], new ApiIntergroup());
							fTarget1 = new MayerGeneralSpherical((P2LennardJones)Species1Potentials[0]);
					        eTarget1 = new MayerESpherical((P2LennardJones)Species1Potentials[0]);
							//fTarget = new MayerGeneral(pIntra1TargetGroup);
			    		    //eTarget = new MayerEGeneral(pIntra1TargetGroup);
						}
						
						
						
	    		 }
	    		 
	    	 }
	    	else{
	    		pInterTargetGroup = new PotentialGroup(2);
	    		if(Species1IntraNonBondedPotentialFlag){
	    			pIntra1TargetGroup = new PotentialGroup(2);
	    		}
	    		if(Species2IntraNonBondedPotentialFlag){
	    			pIntra2TargetGroup = new PotentialGroup(2);
	    		}
	    	}
	    	 
	    	 
	    	 
	    	 
	    	 
	    } catch (SecurityException e) {
	    	 // TODO Auto-generated catch block
	    	 e.printStackTrace();
	     } catch (IllegalArgumentException e) {
	    	 // TODO Auto-generated catch block
	    	 e.printStackTrace();
	     } catch (InstantiationException e) {
	    	 // TODO Auto-generated catch block
	    	 e.printStackTrace();
	     } catch (IllegalAccessException e) {
	    	 // TODO Auto-generated catch block
	    	 e.printStackTrace();
	     } catch (InvocationTargetException e) {
	    	 // TODO Auto-generated catch block
	    	 e.printStackTrace();
	     }
	     
	     space = potential[0].getSpace();
	     ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget1, nPoints>3, eTarget, true);
	     targetCluster.setTemperature(temperature);
	     ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
	     refCluster.setTemperature(temperature);
	     
	     System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
	     
	     final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,potential[0].createSpeciesFactory(), temperature,refCluster,targetCluster);
	     
	     sim.integratorOS.setNumSubSteps(1000);
	     int blocksize = 100;
	     sim.setAccumulatorBlockSize(blocksize);
	        
	       
	        
	        
	        if (false) {
	            double size = 5;
	            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
	            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
	            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
	            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
	            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
	            displayBox0.setPixelUnit(new Pixel(300.0/size));
	            displayBox1.setPixelUnit(new Pixel(300.0/size));
	            displayBox0.setShowBoundary(false);
	            displayBox1.setShowBoundary(false);
	            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
	            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);
	            
	            
	            DiameterHashByType diameterManager = (DiameterHashByType)displayBox0.getDiameterHash();
	           // diameterManager.setDiameter(typeCH2, 0.8*sigmaCH2);
	            diameterManager.setDiameter(sim.getSpecies(0).getAtomType(0),1.0);
	            displayBox1.setDiameterHash(diameterManager);
	            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
	            displayBox0.setColorScheme(colorScheme);
	            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
	            displayBox1.setColorScheme(colorScheme);
	            simGraphic.makeAndDisplayFrame();

	            sim.integratorOS.setNumSubSteps(1000);
	            sim.setAccumulatorBlockSize(1000);
	                
	            // if running interactively, set filename to null so that it doens't read
	            // (or write) to a refpref file
	            sim.getController().removeAction(sim.ai);
	            sim.getController().addAction(new IAction() {
	                public void actionPerformed() {
	                    sim.initRefPref(null, 10);
	                    sim.equilibrate(null, 20);
	                    sim.ai.setMaxSteps(Long.MAX_VALUE);
	                }
	            });
	            sim.getController().addAction(sim.ai);
	            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
	                throw new RuntimeException("Oops");
	            }
	            
	            final DisplayTextBox averageBox = new DisplayTextBox();
	            averageBox.setLabel("Average");
	            final DisplayTextBox errorBox = new DisplayTextBox();
	            errorBox.setLabel("Error");
	            JLabel jLabelPanelParentGroup = new JLabel("B"+nPoints+" (L/mol)^"+(nPoints-1));
	            final JPanel panelParentGroup = new JPanel(new java.awt.BorderLayout());
	            panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
	            panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.WEST);
	            panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
	            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());
	            
	            IAction pushAnswer = new IAction() {
	                public void actionPerformed() {
	                    double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
	                    double ratio = ratioAndError[0];
	                    double error = ratioAndError[1];
	                    data.x = ratio;
	                    averageBox.putData(data);
	                    data.x = error;
	                    errorBox.putData(data);
	                }
	                
	                DataDouble data = new DataDouble();
	            };
	            IEtomicaDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
	            Unit unit = new CompoundUnit(new Unit[]{new UnitRatio(Liter.UNIT, Mole.UNIT)}, new double[]{nPoints-1});
	            averageBox.putDataInfo(dataInfo);
	            averageBox.setLabel("average");
	            averageBox.setUnit(unit);
	            errorBox.putDataInfo(dataInfo);
	            errorBox.setLabel("error");
	            errorBox.setPrecision(2);
	            errorBox.setUnit(unit);
	            sim.integratorOS.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));
	            
	            return;
	        }
	        
	        IAction progressReport = new IAction() {
	            public void actionPerformed() {
	                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
	                double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
	               // System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
	            }
	        };
	        IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
	        progressReportListener.setInterval((int)(steps/10));
	        sim.integratorOS.getEventManager().addListener(progressReportListener);

	        sim.ai.setMaxSteps(steps);
	        for (int i=0; i<2; i++) {
	            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize());
	        }
	        sim.getController().actionPerformed();

	        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
	        System.out.println("actual reference step frequency "+sim.integratorOS.getActualStepFreq0());
	        
	        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
	        System.out.println("ratio average: "+ratioAndError[0]+", error: "+ratioAndError[1]);
	        //System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
	        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(0);
	        System.out.println("hard sphere ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
	                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
	        System.out.println("hard sphere   average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
	                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
	                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
	        System.out.println("hard sphere overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
	                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
	                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);

	        System.out.println("hard sphere autocorrelation function: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.BLOCK_CORRELATION.index)).getData()[0]);
	        
	        System.out.println("hard sphere overlap autocorrelation function: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.BLOCK_CORRELATION.index)).getData()[1]);
	        
	        System.out.println();
	        
	        
	        allYourBase = (DataGroup)sim.accumulators[1].getData(0);
	        System.out.println("lennard jones ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO.index)).getData()[1]
	                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverage.StatType.RATIO_ERROR.index)).getData()[1]);
	        System.out.println("lennard jones average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
	                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
	                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
	        System.out.println("lennard jones overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
	                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
	                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);

	        System.out.println("lennard jones autocorrelation function: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.BLOCK_CORRELATION.index)).getData()[0]);
	        
	        System.out.println("lennard jones overlap autocorrelation function: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.BLOCK_CORRELATION.index)).getData()[1]);
	    
	    
		
	}
	public int Factorial(int n)
	{
		if (n == 0)
			return 1;
		else
			return n * Factorial(n-1);
	}

	

}
