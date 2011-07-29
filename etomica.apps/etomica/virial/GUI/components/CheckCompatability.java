package etomica.virial.GUI.components;

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
	
	
	@SuppressWarnings("unchecked")
	public void checkIfCompatible(ParameterMapping potential1,ParameterMapping potential2,SimulationEnvironment SimENV){
		SimEnv = SimENV;
		potential = new ParameterMapping[2];
		potential[0] = potential1;
		potential[1] = potential2;

		
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
						Species1IntraBondedPotentialFlag = false;
						pIntra1TargetGroup = null;
					}
					if(i==1){
						Species2AtomicFlag = false;
						Species2IntraBondedPotentialFlag = false;
						pIntra2TargetGroup = null;
					}
				}
				if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom(potential[i].getPotential())){
					
						//Inter-bonded potentials is to be gathered especially for alkanes, alcohol, etc
						if(i==0){
							Species1AtomicFlag = true;
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
			
			if(Species1AtomicFlag && Species2AtomicFlag){
				if(potential1.getPotential().equals(potential2.getPotential())){
					potentialCheckFlag = true;
					SameNonBondedPotentialFlag = true;
					
					
				}
				else{
					potentialCheckFlag = false;
				}
			}
			else if(!Species1AtomicFlag && Species2AtomicFlag){
				if(potential1.getPotential().equals(potential2.getPotential())){
					potentialCheckFlag = true;
					
				}
				else if( potential1.getNonBondedInteractionModel().equals(potential2.getNonBondedInteractionModel())){
					potentialCheckFlag = true;
				}
				else{
					potentialCheckFlag = false;
				}
			}
			else if(Species1AtomicFlag && !Species2AtomicFlag){
				if(potential1.getPotential().equals(potential2.getPotential())){
					potentialCheckFlag = true;
				}
				else if( potential1.getNonBondedInteractionModel().equals(potential2.getNonBondedInteractionModel()) ){
					potentialCheckFlag = true;
				}
				else{
					potentialCheckFlag = false;
				}
			}
			else if(!Species1AtomicFlag && !Species2AtomicFlag){
				if(potential1.getPotential().equals(potential2.getPotential())){
					potentialCheckFlag = true;	
					
					if(potential1.getPotential().getConstructors().length > 1){
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
							for(int j=0;j<potential[0].getParametersArray().length;j++){
								ParamValueObj[j+1]=potential[i].getDoubleDefaultParameters(potential[i].getParametersArray()[j].toUpperCase()+potential[i].getPotentialSiteAtIndex(0));
								ParamValueCrossObj[0][j]=potential[i].getParametersArray()[j].toUpperCase();
								ParamValueCrossObj[i+1][j]=potential[i].getDoubleDefaultParameters(potential[i].getParametersArray()[j].toUpperCase()+potential[i].getPotentialSiteAtIndex(0));
							}
							
							
							try{
								if(i==0){
									Species1Molecular = potential[i].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
									
								}
								else if(i==1){
									Species2Molecular = potential[i].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
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
							CrossPotential = potential1.getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
						} catch (Exception e) {
							
							e.printStackTrace();
						}
						
					}
					else if(potential1.getPotential().getConstructors().length == 1){
						
					}
					
				}
				else if( potential1.getNonBondedInteractionModel().equals(potential2.getNonBondedInteractionModel()) ){
					potentialCheckFlag = true;
				}
				else{
					potentialCheckFlag = false;
				}
			}
			
			

			if(potentialCheckFlag){
				//Condition for electrostatic interaction included
				if((Species1ChargeFlag && Species2ChargeFlag)){
					//run simulation
					System.out.println("Run");
					electrostaticPotentialCheckFlag = true;
				}
				else if(Species1MomentFlag && Species2MomentFlag){
					//run simulation
					System.out.println("Run");
					electrostaticPotentialCheckFlag = true;
				}
				else if(!Species1MomentFlag && !Species2MomentFlag && !Species1ChargeFlag && !Species2ChargeFlag){
					//run simulation
					System.out.println("Run");
					electrostaticPotentialCheckFlag = true;
				}
				//-----------------No charge or moment on species 1 but charge or moment on species 2
				
				else if(!Species1ChargeFlag && !Species1MomentFlag && Species2MomentFlag){
					//run simulation
					System.out.println("Run");
					electrostaticPotentialCheckFlag = true;
				}
				
				else if(!Species1ChargeFlag && !Species1MomentFlag && Species2ChargeFlag ){
					//run simulation
					System.out.println("Run");
					electrostaticPotentialCheckFlag = true;
				}
				
				//------------------Charge or moment on species 1 but no charge or moment on species 2
				else if(Species1MomentFlag && !Species2ChargeFlag && !Species2MomentFlag){
					//run simulation
					System.out.println("Run");
					electrostaticPotentialCheckFlag = true;
				}
				
				else if(Species1ChargeFlag & !Species2ChargeFlag && !Species2MomentFlag){
					//run simulation
					System.out.println("Run");
					electrostaticPotentialCheckFlag = true;
				}
				//--------------------------------Charge on species 1 and Moment on species 2 || Moment on species 1 and charge on species 2
				else if(Species1ChargeFlag && Species2MomentFlag){
					//stop simulation
					System.out.println("Stop");
					electrostaticPotentialCheckFlag = false;
				}
				else if(Species1MomentFlag && Species2ChargeFlag){

					//stop simulation
					System.out.println("Stop");
					electrostaticPotentialCheckFlag = false;
				}
			}
			else{
				//stop simulation
				System.out.println("Stop");
				electrostaticPotentialCheckFlag = false;
			}
			
			

		//end of if statement for potential2 not equal to null	
		}
		
		
	}
}
