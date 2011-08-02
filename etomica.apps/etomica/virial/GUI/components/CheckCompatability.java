package etomica.virial.GUI.components;

import java.lang.reflect.Array;
import java.lang.reflect.InvocationTargetException;

import etomica.api.IPotentialMolecular;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.virial.GUI.containers.DialogBoxPanel;

public class CheckCompatability {

	private DialogBoxPanel messageAlert;
	
	private boolean potentialCheckFlag = false;
	private boolean electrostaticPotentialCheckFlag = false;
	private boolean SameNonBondedPotentialFlag = false;
	
	private boolean Species1AtomicFlag;
	private boolean Species2AtomicFlag;
	
	private int[] ElectrostaticFlag;
	private char[] PotentialTypeFlag;
	private boolean Species1ChargeFlag = false;
	private boolean Species2ChargeFlag = false;
	
	private boolean Species1MomentFlag = false;
	private boolean Species2MomentFlag = false;
	
	private boolean InterNonBondedPotentialFlag;
	private boolean Species1IntraBondedPotentialFlag;
	private boolean Species2IntraBondedPotentialFlag;
	
	private PotentialGroup pIntra1TargetGroup;
	private PotentialGroup pIntra2TargetGroup;
	
	
	private SimulationEnvironment SimEnv;
	private ParameterMapping[] potential;
	
	private Object Species1Molecular;
	private Object Species2Molecular;
	private Object CrossPotential;
	
	private Object[] AllPotentials;
	
	@SuppressWarnings("unchecked")
	public void checkIfCompatible(ParameterMapping potential1,ParameterMapping potential2,SimulationEnvironment SimENV){
		SimEnv = SimENV;
		potential = new ParameterMapping[2];
		potential[0] = potential1;
		potential[1] = potential2;

		PotentialTypeFlag = new char[2];
		
		//Are we having a mixture? 
		if(potential2 != null){
			//Yes, We have a mixture of species!!!
			InterNonBondedPotentialFlag = true;
		}
		else{
			InterNonBondedPotentialFlag = false;
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
						Species1AtomicFlag = false;
						PotentialTypeFlag[0] = 'M';
						Species1IntraBondedPotentialFlag = false;
						pIntra1TargetGroup = null;
					}
					if(i==1){
						Species2AtomicFlag = false;
						PotentialTypeFlag[1] = 'M';
						Species2IntraBondedPotentialFlag = false;
						pIntra2TargetGroup = null;
					}
				}
				if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom(potential[i].getPotential())){
					
						//Inter-bonded potentials is to be gathered especially for alkanes, alcohol, etc
						if(i==0){
							Species1AtomicFlag = true;
							PotentialTypeFlag[0] = 'A';
							if(SimEnv.getAlkane1Spheres() > 1){
								Species1IntraBondedPotentialFlag = true;
							}
							else{
								Species1IntraBondedPotentialFlag = false;
								pIntra1TargetGroup = null;
							}
						}
						if(i==1){
							Species2AtomicFlag = true;
							PotentialTypeFlag[1] = 'A';
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
			
			
			if(PotentialTypeFlag[0] == PotentialTypeFlag[1]){
				//Both are molecular Potentials or Both are Atomic Potentials
				if(PotentialTypeFlag[0] == 'M'){
					if(potential1.getPotential().equals(potential2.getPotential())){
						potentialCheckFlag = true;	
						makeMolecularPairPotentials(potential,1);
						
					}
					else if( potential1.getNonBondedInteractionModel().equals(potential2.getNonBondedInteractionModel()) ){
						potentialCheckFlag = true;
						makeMolecularPairPotentials(potential,2);
					}
				}
				
			}
			else{
				//If Species1 is molecular Potential and Species 2 is Atomic and vice versa
				if( potential1.getNonBondedInteractionModel().equals(potential2.getNonBondedInteractionModel()) ){
					potentialCheckFlag = true;
				}
				else{
					potentialCheckFlag = false;
				}
			}
				
			if(potentialCheckFlag){
				//Condition for electrostatic interaction included
				ElectrostaticFlag = getElectrostatic();
				if(ElectrostaticFlag[0] == ElectrostaticFlag[1] || ElectrostaticFlag[0] == 0 || ElectrostaticFlag[1] == 0){
					System.out.println("Run");
				}
				else{
					//stop simulation
					System.out.println("Stop");
				}
				
			}
			else{
				//stop simulation
				System.out.println("Stop");
				
			}
			
			

		//end of if statement for potential2 not equal to null	
		}
		
		
	}
	
	
	
	@SuppressWarnings("unchecked")
	public void makeMolecularPairPotentials(ParameterMapping[] potential, int type){
		boolean OneConstructor = false;
		if(type == 1){
			
			ParameterMapping potential1 = potential[0];
			if(potential1.getPotential().getConstructors().length > 1){
				OneConstructor = false;
			}
			else if(potential1.getPotential().getConstructors().length == 1){
				OneConstructor = true;
			}
			int numberofParameters = potential1.getParametersArray().length;
			Object[][] ParamValueCrossObj = new Object[numberofParameters][3];
			
			@SuppressWarnings("rawtypes")
			Class[] ParamClass = new Class[numberofParameters+1];
			ParamClass[0] = ISpace.class;
			for(int j = 1;j<=numberofParameters;j++){
				ParamClass[j] = Double.TYPE;
			}
			
			for(int i = 0;i < 2;i++){					
				Object[] ParamValueObj = new Object[numberofParameters+1];

				ParamValueObj[0] = potential[0].getSpace();
				if(!OneConstructor){
					for(int j=0;j<potential[0].getParametersArray().length;j++){
						ParamValueObj[j+1]=potential[i].getDoubleDefaultParameters(potential[i].getParametersArray()[j].toUpperCase()+potential[i].getPotentialSiteAtIndex(0));
						ParamValueCrossObj[0][j]=potential[i].getParametersArray()[j].toUpperCase();
						ParamValueCrossObj[i+1][j]=potential[i].getDoubleDefaultParameters(potential[i].getParametersArray()[j].toUpperCase()+potential[i].getPotentialSiteAtIndex(0));
					}
				}
				
				
				try{
					if(i==0){
						if(OneConstructor){
							Species1Molecular = potential[i].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
						}else{
							Species1Molecular = potential[i].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
						}
						
						
					}
					else if(i==1){
						if(OneConstructor){
							Species2Molecular = potential[i].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
						}
						else{
							Species2Molecular = potential[i].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
						}
						
							
					}
				}
				catch(Exception E){
					E.printStackTrace();
				}
			}
			Object[] ParamValueObj = new Object[numberofParameters+1];

			ParamValueObj[0] = potential[0].getSpace();
			String valueA;
			String valueB;
			for(int i=0;i < ParamValueCrossObj.length;i++){
				valueA = ParamValueCrossObj[1][i].toString();
				double ValueA = Double.parseDouble(valueA);
				valueB = ParamValueCrossObj[2][i].toString();
				double ValueB = Double.parseDouble(valueB);
				if(ParamValueCrossObj[0][i].toString().contains("SIGMA")){
					ParamValueObj[i+1] = 0.5*(ValueA+ValueB);
				}
				if(ParamValueCrossObj[0][i].toString().contains("EPSILON")){
					ParamValueObj[i+1] = Math.sqrt(ValueA*ValueB);
				}
				if(ParamValueCrossObj[0][i].toString().contains("MOMENT")){
					ParamValueObj[i+1] = ValueA*ValueB;
				}
				if(ParamValueCrossObj[0][i].toString().contains("MOMENTSQR")){
					ParamValueObj[i+1] = ValueA*ValueB;
				}
			}
			try {
				if(OneConstructor){
					CrossPotential = potential1.getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
				}else{
					CrossPotential = potential1.getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
				}
				
			} catch (Exception e) {
				e.printStackTrace();
			}
			
		}
		if(type == 2){
			ParameterMapping potential1 = potential[0];
			ParameterMapping potential2 = potential[1];
			if(potential1.getNonBondedInteractionModel() == "LennardJones"){
				String[] P1 = potential1.getPotentialSites();
				String[] P2 = potential2.getPotentialSites();
				for(int k = 0;k < 2;k++){
					if(potential[k].getPotential().getConstructors().length > 1){
						OneConstructor = false;
					}
					else if(potential[k].getPotential().getConstructors().length == 1){
						OneConstructor = true;
					}
					int numberofParameters = potential[k].getParametersArray().length;
					Object[][] ParamValueCrossObj = new Object[numberofParameters][3];
				
					@SuppressWarnings("rawtypes")
					Class[] ParamClass = new Class[numberofParameters+1];
					ParamClass[0] = ISpace.class;
					for(int j = 1;j<=numberofParameters;j++){
						ParamClass[j] = Double.TYPE;
					}
				
					Object[] ParamValueObj = new Object[numberofParameters+1];
					ParamValueObj[0] = potential[0].getSpace();
					if(!OneConstructor){
						for(int j=0;j<potential[0].getParametersArray().length;j++){
							ParamValueObj[j+1]=potential[k].getDoubleDefaultParameters(potential[k].getParametersArray()[j].toUpperCase()+potential[k].getPotentialSiteAtIndex(0));
							ParamValueCrossObj[0][j]=potential[k].getParametersArray()[j].toUpperCase();
							ParamValueCrossObj[k+1][j]=potential[k].getDoubleDefaultParameters(potential[k].getParametersArray()[j].toUpperCase()+potential[k].getPotentialSiteAtIndex(0));
						}
					}
					try{
						if(k==0){
							if(OneConstructor){
								Species1Molecular = potential[k].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
							}else{
								Species1Molecular = potential[k].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
							}
						}
						else if(k==1){
							if(OneConstructor){
								Species2Molecular = potential[k].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
							}
							else{
								Species2Molecular = potential[k].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
							}	
						}
					}
					catch(Exception E){
						E.printStackTrace();
					}
					
				}
				
				
			}
			if(potential1.getNonBondedInteractionModel() == "Exp-6"){
				
			}
			
		}
	
	}
	/*
	 * Function to find the charges/moments on the potentials
	 * Input: None
	 * Output: integer Array
	 */
	public int[] getElectrostatic(){
		if((Species1ChargeFlag && Species2ChargeFlag)){
			//run simulation
			
			return new int[]{1,1};
		}
		else if(Species1MomentFlag && Species2MomentFlag){
			//makeMolecularPairPotentials(potential);
			return new int[]{2,2};
		}
		else if(!Species1MomentFlag && !Species2MomentFlag && !Species1ChargeFlag && !Species2ChargeFlag){
			return new int[]{0,0};
		}
		//-----------------No charge or moment on species 1 but charge or moment on species 2
		
		else if(!Species1ChargeFlag && !Species1MomentFlag && Species2MomentFlag){
			return new int[]{0,2};
		}
		
		else if(!Species1ChargeFlag && !Species1MomentFlag && Species2ChargeFlag ){
			return new int[]{0,1};
		}
		
		//------------------Charge or moment on species 1 but no charge or moment on species 2
		else if(Species1MomentFlag && !Species2ChargeFlag && !Species2MomentFlag){
			return new int[]{2,0};
		}
		
		else if(Species1ChargeFlag & !Species2ChargeFlag && !Species2MomentFlag){
			return new int[]{1,0};
		}
		//--------------------------------Charge on species 1 and Moment on species 2 || Moment on species 1 and charge on species 2
		else if(Species1ChargeFlag && Species2MomentFlag){
			return new int[]{1,2};
		}
		else if(Species1MomentFlag && Species2ChargeFlag){
			return new int[]{2,1};
		}
		else{
		return null;
		}
	}
	
	
	
}
