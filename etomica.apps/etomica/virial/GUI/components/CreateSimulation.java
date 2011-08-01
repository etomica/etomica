package etomica.virial.GUI.components;

import java.awt.Color;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import etomica.action.IAction;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiIntergroup;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2Spherical;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.units.CompoundDimension;
import etomica.units.CompoundUnit;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;
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
import etomica.virial.cluster.Standard;
import etomica.virial.simulations.SimulationVirialOverlap;

public class CreateSimulation {
	
	
	
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
	private boolean MolecularFlag;
	
	private boolean Species1ChargeFlag;
	private boolean Species2ChargeFlag;
	
	private boolean Species1MomentFlag;
	private boolean Species2MomentFlag;
	
	private double temperature;
	private double SigmaHSRef;
	private int steps;
	
	private PotentialGroup p11InterTargetGroup;
	private PotentialGroup p22InterTargetGroup;
	private PotentialGroup p12InterTargetGroup;
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
		/*
		SimEnv = simenv;
		//All Environment variables set first
		temperature = SimEnv.getTemperature();
		steps = SimEnv.getNoOfSteps();
		System.out.println(potential[0].getClass().getName());

		
		//Are we having a mixture? 
		if(potential[1] != null){
			//Yes, We have a mixture of species!!!
			InterNonBondedPotentialFlag = true;
		}
		else{
			InterNonBondedPotentialFlag = false;
			Species2IntraNonBondedPotentialFlag = false;
			Species2IntraBondedPotentialFlag = false;
		}
		
		if(InterNonBondedPotentialFlag == true){
			for(int i = 0;i < 2;i++){
				//We figure details abt each of the potential
				
				//If molecular or atomic
				if (etomica.api.IPotentialMolecular.class.isAssignableFrom(potential[i].getPotential())){

					//If potentialsites existing are greater than 1, although we dont have a intrabonded or intra non-bonded potential, but we have 
					//cross potentials
					if(i==0){
						Species1IntraBondedPotentialFlag = false;
						pIntra1TargetGroup = null;
					}
					if(i==1){
						Species2IntraBondedPotentialFlag = false;
						pIntra2TargetGroup = null;
					}
				}
				if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom(potential[i].getPotential())){
					
						//Inter-bonded potentials is to be gathered especially for alkanes, alcohol, etc
						if(i==0){
							if(SimEnv.getAlkane1Spheres() > 1){
								Species1IntraBondedPotentialFlag = true;
							}
							else{
								Species1IntraBondedPotentialFlag = false;
								pIntra1TargetGroup = null;
							}
						}
						if(i==1){
							if(SimEnv.getAlkane2Spheres() > 1){
								Species2IntraBondedPotentialFlag = true;
							}
							else{
								Species2IntraBondedPotentialFlag = false;
								pIntra2TargetGroup = null;
							}
						}
				}
				
				
				//If potential describing the interaction incluses charges or moments
				String[] tempArray = potential[i].getParametersArray();
				for(int j=0; j< tempArray.length;j++){
					if(tempArray[j].equals("CHARGE")){
						if(i == 0){
							Species1ChargeFlag = true;}
						if(i == 1){
							Species2ChargeFlag = true;}
						
					}
					if(tempArray[j].equals("MOMENT")||tempArray[j].equals("MOMENTSQR")){
						//Dipole /Quadrapole moment exists!
						if(i == 0){
							Species1MomentFlag = true;}
						if(i == 1){
							Species2MomentFlag = true;}

					}
				}
			
				
				
				
			//end of for loop for potentials
			}
			
			
			//Condition for electrostatic interaction included
			if((Species1ChargeFlag && Species2ChargeFlag)|| (Species1MomentFlag && Species2MomentFlag)){
				
				
			}
			
		//end of if statement for potential2 not equal to null	
		}
		
		
		
		
		

		
		
		
		
		
		
		//SigmaHSRef will vary according to mixing rules
		SigmaHSRef = SimEnv.getSigmaHSRef();
		if(!InterNonBondedPotentialFlag){
			
		}
		//For Alkane mixtures
		
		final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(SigmaHSRef);
        HSB[3] = Standard.B3HS(SigmaHSRef);
        HSB[4] = Standard.B4HS(SigmaHSRef);
        HSB[5] = Standard.B5HS(SigmaHSRef);
        HSB[6] = Standard.B6HS(SigmaHSRef);
        HSB[7] = Standard.B7HS(SigmaHSRef);
        HSB[8] = Standard.B8HS(SigmaHSRef);
		
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
							pIntra1TargetGroup = new PotentialGroup(2);
							pIntra1TargetGroup.addPotential((IPotentialAtomic)Species1Potentials[0], new ApiIntergroup());
							fTarget = new MayerGeneral(pIntra1TargetGroup);
			    		    eTarget = new MayerEGeneral(pIntra1TargetGroup);
						}
						
						
						
	    		 }
	    		 
	    	 }
	    	else{
	    		p11InterTargetGroup = new PotentialGroup(2);
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
	     
	     ClusterAbstract targetCluster = Standard.virialCluster(nPoints, fTarget, nPoints>3, eTarget, true);
	     targetCluster.setTemperature(temperature);
	     ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
	     refCluster.setTemperature(temperature);
	     
	     System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
	     
	     final SimulationVirialOverlap sim = new SimulationVirialOverlap(space,potential[0].createSpeciesFactory(), temperature,refCluster,targetCluster);
	     
	     sim.integratorOS.setNumSubSteps(1000);
	     int blocksize = 100;
	     sim.setAccumulatorBlockSize(blocksize);
	        
	        IAction progressReport = new IAction() {
	            public void actionPerformed() {
	                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
	                double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
	                System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
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

	        
		
	}
	public int Factorial(int n)
	{
		if (n == 0)
			return 1;
		else
			return n * Factorial(n-1);
	}

	*/
	}

}
