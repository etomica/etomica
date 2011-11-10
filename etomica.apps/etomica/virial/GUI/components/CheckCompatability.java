package etomica.virial.GUI.components;

import java.awt.List;
import java.lang.reflect.Array;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;

import etomica.api.IAtomType;
import etomica.api.IPotentialMolecular;
import etomica.api.ISpecies;
import etomica.atom.iterator.ApiBuilder;
import etomica.potential.P2Exp6Buckingham;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.virial.MCMoveClusterTorsionMulti;
import etomica.virial.GUI.containers.DialogBoxPanel;
import etomica.virial.GUI.models.Sites;

public class CheckCompatability {

	private DialogBoxPanel messageAlert;
	
	//Main Flag to continue Simulation after potentials match
	private boolean potentialCheckFlag = false;
	
	//Flag for Electrostatic Potentials
	private boolean electrostaticPotentialCheckFlag = false;
	
	
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

	
	private SimulationEnvironment SimEnv;
	private ParameterMapping[] potential;
	
	
	
	
	private ArrayList<String> SiteNotParticipatingIn1 = new ArrayList<String>();
	private ArrayList<String[]> UnlikePairsWithoutDuplicates;
	private ArrayList<String[]> UnlikePairsWithoutCrossPairs;
	
	private PotentialObject PObject;
	

	@SuppressWarnings("unchecked")
	public PotentialObject checkIfCompatible(ParameterMapping potential1,ParameterMapping potential2,SimulationEnvironment SimENV){
		
		PObject = new PotentialObject();
		SimEnv = SimENV;
		potential = new ParameterMapping[2];
		potential[0] = potential1;
		potential[1] = potential2;

		PotentialTypeFlag = new char[2];
		
		//Are we having a mixture? 
		if(potential2 != null){
			//Yes, We have a mixture of species!!!
			InterNonBondedPotentialFlag = true;
			PObject.setpInterGroupIJ(null);
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
						
						setSpecies1AtomicFlag(false);
						Species1IntraBondedPotentialFlag = false;
						
					}
					if(i==1){
						setSpecies2AtomicFlag(false);
						Species2IntraBondedPotentialFlag = false;
					}
					PotentialTypeFlag[i] = 'M';
					PObject.setpIntraGroupII(null, i);
					//pIntraTargetGroup[i] = null;
				}
				if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom(potential[i].getPotential())){
					
						//Inter-bonded potentials is to be gathered especially for alkanes, alcohol, etc
						if(i==0){
							setSpecies1AtomicFlag(true);
							PotentialTypeFlag[0] = 'A';
							if(SimEnv.getAlkane1Spheres() > 2){
								Species1IntraBondedPotentialFlag = true;
								PObject.setBondedIIPotentialSetsSpeciesI(new HashMap(), 0);
								PObject.setIteratorSetsSpeciesI(new HashMap(), 0);
								BondedInteraction bInteractions = new BondedInteractionAlkanes(SimEnv.getAlkane1Spheres());
								bInteractions.AddBondedPotentialSets(PObject,potential[i].getSpace(),1);
								if(SimEnv.getAlkane1Spheres() > 3){
									PObject.setTorsionMovesI(new MCMoveClusterTorsionMulti[2]);
								}
							}
							else if(potential[i].getClass().getName().contains("Methanol") || potential[i].getClass().getName().contains("Ethanol")){
								PObject.setBondedIIPotentialSetsSpeciesI(new HashMap(), 0);
								PObject.setIteratorSetsSpeciesI(new HashMap(), 0);
								Species1IntraBondedPotentialFlag = true;
								BondedInteraction bInteractions = new BondedInteractionAlcohols();
								bInteractions.AddBondedPotentialSets(PObject,potential[i].getSpace(),1);
							}
							else{
								Species1IntraBondedPotentialFlag = false;
								PObject.setpIntraGroupII(null, i);
								//pIntraTargetGroup = null;
							}
						}
						if(i==1){
							setSpecies1AtomicFlag(true);
							PotentialTypeFlag[1] = 'A';
							if(SimEnv.getAlkane2Spheres() > 2){
								
								Species2IntraBondedPotentialFlag = true;
								PObject.setBondedIIPotentialSetsSpeciesI(new HashMap(), 1);
								PObject.setIteratorSetsSpeciesI(new HashMap(), 1);
								BondedInteraction bInteractions = new BondedInteractionAlkanes(SimEnv.getAlkane2Spheres());
								bInteractions.AddBondedPotentialSets(PObject,potential[i].getSpace(),2);
								if(SimEnv.getAlkane2Spheres() > 3){
									PObject.setTorsionMovesJ(new MCMoveClusterTorsionMulti[2]);
								}
							}else if(potential[i].getClass().getName().contains("Methanol") || potential[i].getClass().getName().contains("Ethanol")){
								
								Species2IntraBondedPotentialFlag = true;
								PObject.setBondedIIPotentialSetsSpeciesI(new HashMap(), 1);
								PObject.setIteratorSetsSpeciesI(new HashMap(), 1);
								BondedInteraction bInteractions = new BondedInteractionAlcohols();
								bInteractions.AddBondedPotentialSets(PObject,potential[i].getSpace(),2);
							}
							else{
								Species2IntraBondedPotentialFlag = false;
								PObject.setpIntraGroupII(null, i);
								//pIntra2TargetGroup = null;
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
						//The cross potential is also of same type
						makeMolecularPairPotentials(potential,1);
						
					}
					else if( potential1.getNonBondedInteractionModel().equals(potential2.getNonBondedInteractionModel()) ){
						//CrossPotentials are site-site interactions and modelled using either lennard Jones or Exp-6 potential
						potentialCheckFlag = true;
						makeMolecularPairPotentials(potential,2);
					}
				}
				
				
				if(PotentialTypeFlag[0] == 'A'){
					
					if(potential1.getPotential().equals(potential2.getPotential())){
						potentialCheckFlag = true;
						makeAtomicPairPotentials(potential,1);

					}
					else if( potential1.getNonBondedInteractionModel().equals(potential2.getNonBondedInteractionModel()) ){
						potentialCheckFlag = true;
						makeAtomicPairPotentials(potential,1);
					}
				}
			}
			else{
				//If Species1 is molecular Potential and Species 2 is Atomic and vice versa
				if( potential1.getNonBondedInteractionModel().equals(potential2.getNonBondedInteractionModel()) ){
					potentialCheckFlag = true;
					if(PotentialTypeFlag[0] == 'M'){
						ParameterMapping[] TempPotential = new ParameterMapping[2];
						TempPotential[0] = potential[0];
						makeMolecularPairPotentials(TempPotential,3);
						//Will make like pair interactions and unlike pair interactions
						makeAtomicPairPotentials(potential,3);
						
					}
					else{
						ParameterMapping[] TempPotential = new ParameterMapping[2];
						TempPotential[0] = potential[1];
						makeMolecularPairPotentials(TempPotential,3);
						//Will make like pair interactions and unlike pair interactions
						makeAtomicPairPotentials(potential,2);
					}

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
					return PObject;
				}
				else{
					//stop simulation
					System.out.println("Stop");
					return null;
				}
				
			}
			else{
				//stop simulation
				System.out.println("Stop");
				return null;
				
			}
		//end of if statement for potential2 not equal to null	
		}
		else{
			return null;
		}
		
		
		
	}
	
	
	
	private void makeAtomicPairPotentials(ParameterMapping[] potential, int type) {
		// TODO Auto-generated method stub
		ParameterMapping potential1 = potential[0];
		ParameterMapping potential2 = potential[1];
		
		//Get the List of Sites from each of the species
		String[] NewP1 = getSitesNotParticipating(potential1);
		String[] NewP2 = getSitesNotParticipating(potential2);
		
		//Type 1 is where both are atomic potentials // AtomSets are created separately for each set of likepair interactions
		if(type == 1){
			
			ArrayList<String[]> AllPairs = new ArrayList<String[]>();
	
			//For like interactions
			for(int i=0;i<2;i++){
				
				String[] PotentialSites1 = getSitesNotParticipating(potential[i]);
				String[] PotentialSites2 = PotentialSites1;
				
				// Make Pair Sites
				ArrayList<String[]> TempPairSites = makePairSites(PotentialSites1,PotentialSites2);
				
				UnlikePairsWithoutDuplicates = removeDulplicatesFromPairSites(TempPairSites);
				
				//Removes Similar Cross Pairs Example : C,H and H,C
				UnlikePairsWithoutCrossPairs = removeSimilarCrossPairSites(UnlikePairsWithoutDuplicates);
				
				
				ParameterMapping[] PotentialTemp = new ParameterMapping[2];
				PotentialTemp[0] = potential[i];
				//AtomSetPure[i] = getPairAtomTypes(UnlikePairsWithoutDuplicates,PotentialTemp);
				
				PObject.setAtomSetPureSpeciesI(getPairAtomTypes(UnlikePairsWithoutDuplicates,PotentialTemp), i);
				for(int j=0;j<UnlikePairsWithoutDuplicates.size();j++){
					String[] tempString = UnlikePairsWithoutDuplicates.get(j);
					AllPairs.add(tempString);
				}
				
			}
			
			
			// Make Pair Sites for UnLike - Interactions
			ArrayList<String[]> TempPairSites = makePairSites(NewP1,NewP2);
			
			UnlikePairsWithoutDuplicates = removeDulplicatesFromPairSites(TempPairSites);
			
			for(int j=0;j<UnlikePairsWithoutDuplicates.size();j++){
				String[] tempString = UnlikePairsWithoutDuplicates.get(j);
				AllPairs.add(tempString);
			}
			
			//Removes Similar Cross Pairs Example : C,H and H,C
			UnlikePairsWithoutCrossPairs = removeSimilarCrossPairSites(UnlikePairsWithoutDuplicates);
			
			//AtomSetMix = getPairAtomTypes(UnlikePairsWithoutDuplicates,potential);
			PObject.setAtomSetMix(getPairAtomTypes(UnlikePairsWithoutDuplicates,potential));
			
			//Creating a common pool of potentials to avoid repetition!
			UnlikePairsWithoutDuplicates = removeDulplicatesFromPairSites(AllPairs);
			
			UnlikePairsWithoutCrossPairs = removeSimilarCrossPairSites(UnlikePairsWithoutDuplicates);
			
			if(potential[0].getNonBondedInteractionModel() == "LennardJones"){
				//PotentialSets = getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class);
				PObject.setPotentialSets(getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class));
			}
			
		}
		
		//Type 2 and 3 is where one of the potentials is molecular
		if(type == 2 || type == 3){
		
			String[] PotentialSites1 = getSitesNotParticipating(potential[type-2]);
			String[] PotentialSites2 = PotentialSites1;
			
			ArrayList<String[]> AllPairs = new ArrayList<String[]>();
			// Make Pair Sites
			ArrayList<String[]> TempPairSites = makePairSites(PotentialSites1,PotentialSites2);
			
			UnlikePairsWithoutDuplicates = removeDulplicatesFromPairSites(TempPairSites);
			
			//Removes Similar Cross Pairs Example : C,H and H,C
			UnlikePairsWithoutCrossPairs = removeSimilarCrossPairSites(UnlikePairsWithoutDuplicates);
			
			
			ParameterMapping[] PotentialTemp = new ParameterMapping[2];
			PotentialTemp[0] = potential[type-2];
			//AtomSetPure[0] = getPairAtomTypes(UnlikePairsWithoutDuplicates,PotentialTemp);
			PObject.setAtomSetPureSpeciesI(getPairAtomTypes(UnlikePairsWithoutDuplicates,PotentialTemp), 0);
			
			for(int j=0;j<UnlikePairsWithoutDuplicates.size();j++){
				String[] tempString = UnlikePairsWithoutDuplicates.get(j);
				AllPairs.add(tempString);
			}

			getSameSitesOnBothSpecies(NewP1,NewP2,potential);
			
			// Make Unlike Pair Sites
			ArrayList<String[]> CrossPairSites = makePairSites(NewP1,NewP2);
			
			//Remove Duplicates
			UnlikePairsWithoutDuplicates = removeDulplicatesFromPairSites(CrossPairSites);
			
			//Removes Similar Cross Pairs Example : C,H and H,C
			UnlikePairsWithoutCrossPairs = removeSimilarCrossPairSites(UnlikePairsWithoutDuplicates);
			
			//Get the Atom Sets for All Pair Sites including Similar Cross pairs
			//AtomSetMix = getPairAtomTypes(UnlikePairsWithoutDuplicates,potential);
			PObject.setAtomSetMix(getPairAtomTypes(UnlikePairsWithoutDuplicates,potential));
			
			//Get the Potential Sets for Pair Sites filtered to remove similar cross sites
			if(potential1.getNonBondedInteractionModel() == "LennardJones"){
				PObject.setPotentialSets(getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class));
			}
			
		}
		
		
		
		
		
		
	}



	@SuppressWarnings({ "unchecked", "unused" })
	public void makeMolecularPairPotentials(ParameterMapping[] potential, int type){
		boolean OneConstructor = false;
		
		//Both molecular potentials are same
		if(type == 1){
			
			ParameterMapping potential1 = potential[0];
			OneConstructor = hasSingleConstructorForPotential(potential1.getPotential());
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
				if(!OneConstructor && potential[0].getPotentialSites().length == 1){
					for(int j=0;j<potential[0].getParametersArray().length;j++){
						ParamValueObj[j+1]=potential[i].getDoubleDefaultParameters(potential[i].getParametersArray()[j].toUpperCase()+potential[i].getPotentialSiteAtIndex(0));
						ParamValueCrossObj[0][j]=potential[i].getParametersArray()[j].toUpperCase();
						ParamValueCrossObj[i+1][j]=potential[i].getDoubleDefaultParameters(potential[i].getParametersArray()[j].toUpperCase()+potential[i].getPotentialSiteAtIndex(0));
					}
				}
				
				
				try{
					if(i==0){
						if(OneConstructor){
							PObject.setSpecies1Molecular(potential[i].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]));
							//Species1Molecular = potential[i].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
						}else{
							PObject.setSpecies1Molecular(potential[i].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj));
							//Species1Molecular = potential[i].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
						}
						
						
					}
					else if(i==1){
						if(OneConstructor){
							PObject.setSpecies2Molecular(potential[i].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]));
							//Species2Molecular = potential[i].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
						}
						else{
							PObject.setSpecies2Molecular(potential[i].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj));
							//Species2Molecular = potential[i].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
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
				ParamValueObj[i+1]=GetMixingRules(ValueA,ValueB,ParamValueCrossObj[0][i].toString());
				
			}
			try {
				if(OneConstructor){
					PObject.setCrossPotential(potential1.getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]));
					//CrossPotential = potential1.getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
				}else{
					PObject.setCrossPotential(potential1.getPotential().getConstructor(ParamClass).newInstance(ParamValueObj));
					//CrossPotential = potential1.getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
				}
				
			} catch (Exception e) {
				e.printStackTrace();
			}
			
		}
		
		//Molecular potentials are different, but have same non-bonded potential(atomic potential interaction for unlike atoms!!
		if(type == 2){
			ParameterMapping potential1 = potential[0];
			ParameterMapping potential2 = potential[1];
			HashMap PotentialSets = null;
			
				//Molecular potentials for interaction between like atoms
				for(int k = 0;k < 2;k++){
					OneConstructor = hasSingleConstructorForPotential(potential[k].getPotential());
					int numberofParameters = potential[k].getParametersArray().length;
					
				
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
						}
					}
					try{
						if(k==0){
							if(OneConstructor){
								PObject.setSpecies1Molecular(potential[k].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]));
								//Species1Molecular = potential[k].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
								
							}else{
								PObject.setSpecies1Molecular(potential[k].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj));
								//Species1Molecular = potential[k].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
							}
							
						}
						else if(k==1){
							if(OneConstructor){
								PObject.setSpecies2Molecular(potential[k].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]));
								//Species2Molecular = potential[k].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
							}
							else{
								PObject.setSpecies2Molecular(potential[k].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj));
								//Species2Molecular = potential[k].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
							}	
							
						}
					}
					catch(Exception E){
						E.printStackTrace();
					}
					
					
					
				}
				
				//Get the list of sites of each species which do not particpate in Non-Bonded Interactions
				String[] NewP1 = getSitesNotParticipating(potential1);
				String[] NewP2 = getSitesNotParticipating(potential2);

				getSameSitesOnBothSpecies(NewP1,NewP2,potential);
				
				// Make Pair Sites
				ArrayList<String[]> TempPairSites = makePairSites(NewP1,NewP2);
				
				//Remove Duplicates
				UnlikePairsWithoutDuplicates = removeDulplicatesFromPairSites(TempPairSites);
				
				//Removes Similar Cross Pairs Example : C,H and H,C
				UnlikePairsWithoutCrossPairs = removeSimilarCrossPairSites(UnlikePairsWithoutDuplicates);
				
				//Get the Atom Sets for All Pair Sites including Similar Cross pairs
				//HashMap AtomSets = getPairAtomTypes(UnlikePairsWithoutDuplicates,potential);
				PObject.setAtomSetMix(getPairAtomTypes(UnlikePairsWithoutDuplicates,potential));
				
				//Get the Potential Sets for Pair Sites filtered to remove similar cross sites
				
				if(potential1.getNonBondedInteractionModel() == "LennardJones"){
					//PotentialSets = getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class);
					PObject.setPotentialSets(getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class));
				}
				

				if(potential1.getNonBondedInteractionModel() == "Exp-6"){
					//PotentialSets = getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2Exp6Buckingham.class);
					PObject.setPotentialSets(getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2Exp6Buckingham.class));
				}

			}
		
		//Single Molecular Potential to be created
		if(type == 3){
			ParameterMapping potential1 = potential[0];
					//Molecular potentials for interaction between like atoms
					OneConstructor = hasSingleConstructorForPotential(potential[0].getPotential());
					int numberofParameters = potential[0].getParametersArray().length;
					
				
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
							ParamValueObj[j+1]=potential[0].getDoubleDefaultParameters(potential[0].getParametersArray()[j].toUpperCase()+potential[0].getPotentialSiteAtIndex(0));
						}
					}
					try{
							if(OneConstructor){
								//Species1Molecular = potential[0].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
								PObject.setSpecies1Molecular(potential[0].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]));
								
							}else{
								//Species1Molecular = potential[0].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj);
								PObject.setSpecies1Molecular(potential[0].getPotential().getConstructor(ParamClass).newInstance(ParamValueObj));
							}
					}
					catch(Exception E){
						E.printStackTrace();
					}
				
		}

	}
	
	
	/*
	 * Function Retrieves the IAtomTypes for Pair Sites
	 */
	
	
	private void getSameSitesOnBothSpecies(String[] newP1, String[] newP2,
			ParameterMapping[] potential2) {
		
		ArrayList<String> CommonSites= new ArrayList<String>();
		
		
		for(int i=0;i<newP1.length;i++){
			for(int j=0;j<newP2.length;j++){
				if(newP1[i]==newP2[j]){
					if(!getPairSiteValidation(newP1[i],potential2)){
						newP1[i] = newP1[i] + "-a";
						newP2[j] = newP2[j] + "-b";
					}
				}
			}
		}
		
		
		
	}



	@SuppressWarnings({ "unchecked", "rawtypes" })
	private HashMap getPairSitePotentials(
			ArrayList<String[]> unlikePairsWithoutCrossPairs,
			ParameterMapping[] potential,Class TargetPotential) {
		String[] ParametersForPotential = null;
		
		if(TargetPotential.equals(P2LennardJones.class)){
			ParametersForPotential = new String[]{"SIGMA","EPSILON"};
			hasSingleConstructorForPotential(TargetPotential);
			
		}
		
		Double[] Site12 = new Double[ParametersForPotential.length];
		HashMap PairPotentialSets = new HashMap();
		for(int i=0;i<unlikePairsWithoutCrossPairs.size();i++){
			String Site1 = unlikePairsWithoutCrossPairs.get(i)[0];
			String Site2 = unlikePairsWithoutCrossPairs.get(i)[1];
			
			
			
			for(int j=0;j<ParametersForPotential.length;j++){
				
				//if condition to check if SiteA belongs to Species1 selected. If SiteA does not belong to SpeciesA, then SiteA must be part of SpeciesB
				double Value1;
				double Value2;
				
				
				if(potential[1] != null){
					
					
					
					if(Site1 == Site2){
						if(isPotentialSite(Site1,potential[0])){
							if(isPotentialSite(Site1,potential[1])){
								Value1 = potential[0].getDoubleDefaultParameters(ParametersForPotential[j]+Site1);
								Value2 = potential[1].getDoubleDefaultParameters(ParametersForPotential[j]+Site2);
							}
							else{
								Value1 = potential[0].getDoubleDefaultParameters(ParametersForPotential[j]+Site1);
								Value2 = potential[0].getDoubleDefaultParameters(ParametersForPotential[j]+Site2);
							}
						}
						else{
							Value1 = potential[1].getDoubleDefaultParameters(ParametersForPotential[j]+Site1);
							Value2 = potential[1].getDoubleDefaultParameters(ParametersForPotential[j]+Site2);
						}
					}
					else{
						if(Site1.contains("-a") || Site1.contains("-b")){
							String[] SiteA = Site1.split("-");
							Site1 = SiteA[0];
						}
						if(Site2.contains("-a") || Site2.contains("-b")){
							String[] SiteB = Site2.split("-");
							Site2 = SiteB[0];	
						}
						if(isPotentialSite(Site1,potential[0])){
							Value1 = potential[0].getDoubleDefaultParameters(ParametersForPotential[j]+Site1);
							if(isPotentialSite(Site2,potential[0])){
								Value2 = potential[0].getDoubleDefaultParameters(ParametersForPotential[j]+Site2);
							}
							else{
								Value2 = potential[1].getDoubleDefaultParameters(ParametersForPotential[j]+Site2);
							}
						}
						else{
							Value1 = potential[1].getDoubleDefaultParameters(ParametersForPotential[j]+Site1);
							if(isPotentialSite(Site2,potential[0])){
								Value2 = potential[0].getDoubleDefaultParameters(ParametersForPotential[j]+Site2);
							}
							else{
								Value2 = potential[1].getDoubleDefaultParameters(ParametersForPotential[j]+Site2);
							}
						}
						
					}
				}
				else{
					Value1 = potential[0].getDoubleDefaultParameters(ParametersForPotential[j]+Site1);
					Value2 = potential[0].getDoubleDefaultParameters(ParametersForPotential[j]+Site2);
				}
				
				
				
				
					
				Site12[j]=GetMixingRules(Value1,Value2,ParametersForPotential[j]);
			}
			
			Class[] ParamClass = new Class[ParametersForPotential.length+1];
			ParamClass[0] = ISpace.class;
			for(int j = 1;j<=ParametersForPotential.length;j++){
				ParamClass[j] = Double.TYPE;
			}
			
			Object[] ParamValueObj = new Object[ParametersForPotential.length+1];
			ParamValueObj[0] = potential[0].getSpace();
			for(int m=0;m<Site12.length;m++){
				ParamValueObj[m+1]=Site12[m];
			}
			
			try {
				PairPotentialSets.put(unlikePairsWithoutCrossPairs.get(i),TargetPotential.getConstructor(ParamClass).newInstance(ParamValueObj));
			} catch (Exception e) {
				
				e.printStackTrace();
			}
		}
		
		
		
		return PairPotentialSets;
	}


	

	private double GetMixingRules(double value1, double value2, String string) {
		// TODO Auto-generated method stub
		
		double Value = 0.0;
		if(string == "SIGMA"){
			Value = 0.5*(value1+value2);
		}
		
		if(string == "EPSILON"){
			Value =  Math.sqrt(value1+value2);
		}

		if(string == "MOMENT"){
			Value = value1*value2;
		}
		
		if(string == "MOMENTSQR"){
			Value = value1*value2;
		}
		return Value;
		
		
	}



	@SuppressWarnings("unchecked")
	private HashMap getPairAtomTypes(ArrayList<String[]> unLikePairs,
			ParameterMapping[] potential2) {
		ISpecies species1 = potential2[0].createSpeciesFactory();
		Object obj1 = potential2[0].createSpeciesFactory();
		ISpecies species2;
		Object obj2;
		if(potential2[1] != null){
			species2 = potential2[1].createSpeciesFactory();
			 obj2 = potential2[1].createSpeciesFactory();}
		else{
			species2 = species1;
			obj2 = obj1;
		}
		String Site1Name = "";
		String Site2Name = "";
		IAtomType Site1Atom;
		IAtomType Site2Atom;
		HashMap AtomSets = new HashMap();

		for(int i=0;i<unLikePairs.size();i++){
			String Site1 = unLikePairs.get(i)[0];
			String Site2 = unLikePairs.get(i)[1];
			if(Site1.contains("-a") || Site1.contains("-b")){
				String[] SiteA = Site1.split("-");
				Site1 = SiteA[0];
			}
			if(Site2.contains("-a") || Site2.contains("-b")){
				String[] SiteB = Site2.split("-");
				Site2 = SiteB[0];
			}
			Object[] SiteAtoms = new Object[2];
			for(Sites site: Sites.values()){
				if(site.toString().toUpperCase().equals(Site1.toUpperCase())){
					Site1Name = site.getSite();
				}
				if(site.toString().toUpperCase().equals(Site2.toUpperCase())){
					Site2Name = site.getSite();
				}
			}

			Method[] MethodS1 = species1.getClass().getMethods();
			
			for(int j=0;j<MethodS1.length;j++){
				if(MethodS1[j].getReturnType().equals(IAtomType.class)){
					if(MethodS1[j].toString().toUpperCase().contains(Site1.toUpperCase()+"TYPE")|| MethodS1[j].toString().toUpperCase().contains(Site1Name.toUpperCase()+"TYPE")){
						try {
							Site1Atom = (IAtomType) MethodS1[j].invoke(obj1, null);
							SiteAtoms[0] = Site1Atom;
						} catch (Exception e) {
							e.printStackTrace();
						}
						
					}
				}
			}
			
			Method[] MethodS2 = species2.getClass().getMethods();
			
			for(int k=0;k<MethodS2.length;k++){
				if(MethodS2[k].getReturnType().equals(IAtomType.class)){
					if(MethodS2[k].toString().toUpperCase().contains(Site2.toUpperCase()+"TYPE")|| MethodS2[k].toString().toUpperCase().contains(Site2Name.toUpperCase()+"TYPE")){
						System.out.println(Site2+" "+ MethodS2[k].toString());
						try {
							Site2Atom = (IAtomType) MethodS2[k].invoke(obj2, null);
							SiteAtoms[1] = Site2Atom;
						} catch (Exception e) {
							e.printStackTrace();
						} 
					}
				}
			}
			 AtomSets.put(unLikePairs.get(i), SiteAtoms);
		}
		return AtomSets;
		
	}


	
	private ArrayList<String[]> makePairSites(String[] newP1, String[] newP2){
		String[][] PairSitesDraft1on2 = new String[newP1.length*newP2.length][2];
		String[][] PairSitesDraft2on1 = new String[newP2.length*newP1.length][2];
		int index= 0;
		int index1 = 0;
		for(int i=0;i<newP1.length;i++){
			for(int j=0;j<newP2.length;j++){
				PairSitesDraft1on2[index][0]=newP1[i];
				PairSitesDraft1on2[index][1]=newP2[j];
				index++;
			}
		}
		for(int i=0;i<newP2.length;i++){
			for(int j=0;j<newP1.length;j++){
				PairSitesDraft2on1[index1][0]=newP2[i];
				PairSitesDraft2on1[index1][1]=newP1[j];
				index1++;
			}
		}
		
		ArrayList<String[]> PairSitesD = new ArrayList<String[]>();
		for(int i=0;i<PairSitesDraft1on2.length;i++){
			PairSitesD.add(PairSitesDraft1on2[i]);
		}
		for(int i=PairSitesDraft1on2.length;i<PairSitesDraft1on2.length+PairSitesDraft2on1.length;i++){
			PairSitesD.add(PairSitesDraft2on1[i - PairSitesDraft1on2.length]);
		}
		return PairSitesD;
	}
	

	/*
	 * Function to deletes any duplicates 
	 * Input: String Array, String Array, Potentials
	 * Output: ArrayList
	 */
	private ArrayList<String[]> removeDulplicatesFromPairSites(ArrayList<String[]> PairSitesDWithDuplicates) {
		
		
		for(int k=0;k<PairSitesDWithDuplicates.size();k++){
			for(int l=PairSitesDWithDuplicates.size()- 1;l>k;l--){
				String[] P1 = PairSitesDWithDuplicates.get(k);
				String[] P2 = PairSitesDWithDuplicates.get(l);
				if(P1[0]==P2[0]&& P1[1]==P2[1])	{
					/*if(getPairSiteValidation(P1[0],Potential)){
						if(getPairSiteValidation(P1[1],Potential)){
							PairSitesDWithDuplicates.remove(l);
						}
					}*/
					PairSitesDWithDuplicates.remove(l);
				}
				
			}
		}
		return PairSitesDWithDuplicates;
		
	}
	
	
	
	private ArrayList<String[]> removeSimilarCrossPairSites(ArrayList<String[]> PairSitesDWithoutDuplicates) {
		
		ArrayList<String[]> TempArrayList = new ArrayList<String[]>();
		for(int i=0;i<PairSitesDWithoutDuplicates.size();i++){
			TempArrayList.add(PairSitesDWithoutDuplicates.get(i));
		}
		
		for(int k=0;k<TempArrayList.size();k++){
			for(int l=TempArrayList.size()- 1;l>k;l--){
				String[] P1 = TempArrayList.get(k);
				String[] P2 = TempArrayList.get(l);
				if(P1[0]==P2[1]&& P1[1]==P2[0])	{
					/*if(isPotentialSite(P1[0],Potential[0])){
						//Species1 has the site
						
					}else{
						if(isPotentialSite(P1[1],Potential[0]) && isPotentialSite(P1[1],Potential[1])){
							if(getPairSiteValidation(P1[1],Potential)){
								if(getPairSiteValidation(P1[1],Potential)){
									TempArrayList.remove(l);
								}
							}
							
						}
					}
					if(getPairSiteValidation(P1[0],Potential)){
						if(getPairSiteValidation(P1[1],Potential)){
							TempArrayList.remove(l);
						}
					}
					*/
					TempArrayList.remove(l);
				}
				
			}
		}
		return TempArrayList;
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
	
	
	public String[] getSitesNotParticipating(ParameterMapping Potential){
		ArrayList<String> TempArray = new ArrayList<String>();
		String[] NewP;
		for(int i=0;i < Potential.getPotentialSites().length;i++){
			String PotentialSite = Potential.getPotentialSites()[i];
			for(int j=0;j<Potential.getParametersArray().length;j++){
				String Param = Potential.getParametersArray()[j];
				Double Value = Potential.getDoubleDefaultParameters(Param+PotentialSite);
				if(Potential.getNonBondedInteractionModel() =="LennardJones"){
					if(Param == "SIGMA" || Param == "EPSILON"){
						if(Value == 0.0){
							if(!TempArray.contains(PotentialSite)){
								TempArray.add(PotentialSite);
							}
						}
					}
				}
				
				if(Potential.getNonBondedInteractionModel()=="Exp-6"){
					
				}
			}
		}
		if(TempArray.size() > 0){
			NewP = new String[Potential.getPotentialSites().length - TempArray.size()];
			int index = 0;
			for(int i=0;i<Potential.getPotentialSites().length;i++){
				for(int j=0;j<TempArray.size();j++){
					if(Potential.getPotentialSites()[i]!=TempArray.get(j)){
						NewP[index]=Potential.getPotentialSites()[i];
						index++;
					}
				}
			}
			
		}
		else{
			NewP = new String[Potential.getPotentialSites().length];
			for(int m=0;m<Potential.getPotentialSites().length;m++){
				String tempStr = Potential.getPotentialSites()[m].toString();
				NewP[m] = tempStr;
			}
			
		}
		
		return NewP;
		
	}
	
	
	
	
	
	
	public boolean getPairSiteValidation(String PotentialSite, ParameterMapping[] Potential){
		
			Object[][] ParamValueCrossObj = null;
			int ParamCheck = 0;
		
			if(Potential[0].getNonBondedInteractionModel() =="LennardJones"){
				ParamValueCrossObj = new Object[2][3];
				ParamCheck=2;
			}
			
			for(int k=0;k<2;k++){
				for(int l=0;l<Potential[k].getPotentialSites().length;l++){
					if(PotentialSite == Potential[k].getPotentialSiteAtIndex(l)){
						for(int m=0;m<Potential[k].getParametersArray().length;m++){
							String Param = Potential[k].getParametersArray()[m];
							
							if(Param == "SIGMA"|| Param == "EPSILON"){
								ParamValueCrossObj[m][0] = Param;
								if(k==0){
									ParamValueCrossObj[m][1] = Potential[k].getDoubleDefaultParameters(Param+PotentialSite);}
								if(k==1){
									ParamValueCrossObj[m][2] = Potential[k].getDoubleDefaultParameters(Param+PotentialSite);}
							}
						}
					}
				}
			}
			//confirm the values of the common pairs
			int ConfirmIndexFlag=0;
			for(int a=0;a<ParamValueCrossObj.length;a++){
				double ValueP1 = Double.parseDouble(ParamValueCrossObj[a][1].toString());
				double ValueP2 = Double.parseDouble(ParamValueCrossObj[a][2].toString());
				if(ValueP1==ValueP2){
					ConfirmIndexFlag++;
				}
			}
			
			if(ConfirmIndexFlag == ParamCheck){
				return true;
			}
			else{
				return false;
			}
		}
	
	public boolean hasSingleConstructorForPotential(Class PotentialClassName){
		
		Boolean ReturnFlag = false;
		if(PotentialClassName.getConstructors().length>1){
			ReturnFlag = false;
		}else if(PotentialClassName.getConstructors().length==1){
			ReturnFlag = true;
		}
		return ReturnFlag;
		
	}
	
	
	
	public boolean isSpecies1AtomicFlag() {
		return Species1AtomicFlag;
	}



	public void setSpecies1AtomicFlag(boolean species1AtomicFlag) {
		Species1AtomicFlag = species1AtomicFlag;
	}



	public boolean isSpecies2AtomicFlag() {
		return Species2AtomicFlag;
	}



	public void setSpecies2AtomicFlag(boolean species2AtomicFlag) {
		Species2AtomicFlag = species2AtomicFlag;
	}


	public boolean isPotentialSite(String PotentialSite,ParameterMapping Potential){
		
		Boolean ReturnFlag = false;
		for(int i=0;i<Potential.getPotentialSites().length;i++){
			if(PotentialSite==Potential.getPotentialSites()[i]){
				ReturnFlag= true;
				break;
			}
		}
		return ReturnFlag;
		
	}
	
}
