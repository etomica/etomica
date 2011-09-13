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
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.virial.GUI.containers.DialogBoxPanel;
import etomica.virial.GUI.models.Sites;

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
	private PotentialGroup pInterTargetGroup;
	
	
	private SimulationEnvironment SimEnv;
	private ParameterMapping[] potential;
	
	private Object Species1Molecular;
	private Object Species2Molecular;
	private Object CrossPotential;
	private Object[] CrossPotentials;
	
	private Object[][] CrossPotentialsTypes;
	
	private Object[] AllPotentials;
	
	private ArrayList<String> SiteNotParticipatingIn1 = new ArrayList<String>();
	private ArrayList<String[]> UnlikePairsWithoutDuplicates;
	private ArrayList<String[]> UnlikePairsWithoutCrossPairs;
	
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
						
						setSpecies1AtomicFlag(false);
						PotentialTypeFlag[0] = 'M';
						Species1IntraBondedPotentialFlag = false;
						pIntra1TargetGroup = null;
					}
					if(i==1){
						setSpecies2AtomicFlag(false);
						PotentialTypeFlag[1] = 'M';
						Species2IntraBondedPotentialFlag = false;
						pIntra2TargetGroup = null;
					}
				}
				if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom(potential[i].getPotential())){
					
						//Inter-bonded potentials is to be gathered especially for alkanes, alcohol, etc
						if(i==0){
							setSpecies1AtomicFlag(true);
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
							setSpecies1AtomicFlag(true);
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
	
	
	
	@SuppressWarnings({ "unchecked", "unused" })
	public void makeMolecularPairPotentials(ParameterMapping[] potential, int type){
		boolean OneConstructor = false;
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
				ParamValueObj[i+1]=GetMixingRules(ValueA,ValueB,ParamValueCrossObj[0][i].toString());
				
				/*if(ParamValueCrossObj[0][i].toString().contains("SIGMA")){
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
				}*/
				
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
				
				for(int k = 0;k < 2;k++){
					OneConstructor = hasSingleConstructorForPotential(potential[k].getPotential());
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
				
				//Get the list of sites of each species which do not particpate in Non-Bonded Interactions
				String[] NewP1 = getSitesNotParticipating(potential1);
				String[] NewP2 = getSitesNotParticipating(potential2);
				
			
				
				getSameSitesOnBothSpecies(NewP1,NewP2,potential);
				
				
				// Make Pair Sites
				ArrayList<String[]> TempPairSites = makePairSites(NewP1,NewP2);
				
				//Remove Duplicates
				UnlikePairsWithoutDuplicates = removeDulplicatesFromPairSites(TempPairSites, potential);
				
				//Removes Similar Cross Pairs Example : C,H and H,C
				UnlikePairsWithoutCrossPairs = removeSimilarCrossPairSites(UnlikePairsWithoutDuplicates, potential);
				
				//Get the Atom Sets for All Pair Sites including Similar Cross pairs
				HashMap AtomSets = getPairAtomTypes(UnlikePairsWithoutDuplicates,potential);
				
				//Get the Potential Sets for Pair Sites filtered to remove similar cross sites
				HashMap PotentialSets = getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class);
				
				pInterTargetGroup = new PotentialGroup(2);
				
				
				Map AtomSetMap = AtomSets;
				Set AtomEntries = AtomSetMap.entrySet();
				Iterator AtomItr = AtomEntries.iterator();
				
				
				Map PotentialSetsMap = PotentialSets;
				Set PotentialEntries = PotentialSetsMap.entrySet();
				Iterator PotentialItr = PotentialEntries.iterator();
				
				while(PotentialItr.hasNext()){
					Map.Entry PotentialEntry = (Map.Entry) PotentialItr.next();
					String[] PotentialKey= (String[]) PotentialEntry.getKey();
					while (AtomItr.hasNext()) {
						Map.Entry AtomEntry = (Map.Entry) AtomItr.next();
						String[] AtomKey= (String[]) AtomEntry.getKey();
						if(AtomKey[0] == PotentialKey[0] && AtomKey[1] == PotentialKey[1] || AtomKey[0] == PotentialKey[1] && AtomKey[1] == PotentialKey[0]){
							Object[] AtomObjects = (Object[]) AtomEntry.getValue();
							pInterTargetGroup.addPotential((P2LennardJones)PotentialEntry.getValue(),
									ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{(IAtomType)AtomObjects[0],(IAtomType)AtomObjects[1]}));
							System.out.println(pInterTargetGroup.getPotentials().length);
						}
					}
					if(!AtomItr.hasNext()){
						AtomItr = AtomEntries.iterator();
					}
					
                }

				

				
			}
			if(potential1.getNonBondedInteractionModel() == "Exp-6"){
				
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
		//System.out.
		/*if(TargetPotential.equals(P2ElectrostaticWithHardCore.class)){
			ParametersForPotential = new String[]{"CHARGE"};
		}*/
		Double[] Site12 = new Double[ParametersForPotential.length];
		HashMap PairPotentialSets = new HashMap();
		for(int i=0;i<unlikePairsWithoutCrossPairs.size();i++){
			String Site1 = unlikePairsWithoutCrossPairs.get(i)[0];
			String Site2 = unlikePairsWithoutCrossPairs.get(i)[1];
			
			if(Site1.contains("-a") || Site1.contains("-b")){
				String[] SiteA = Site1.split("-");
				Site1 = SiteA[0];
			}
			if(Site2.contains("-a") || Site2.contains("-b")){
				String[] SiteB = Site2.split("-");
				Site2 = SiteB[0];
			}
			
			for(int j=0;j<ParametersForPotential.length;j++){
				double Value1 = potential[0].getDoubleDefaultParameters(ParametersForPotential[j]+Site1);
				double Value2 = potential[1].getDoubleDefaultParameters(ParametersForPotential[j]+Site2);
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
		ISpecies species2 = potential2[1].createSpeciesFactory();
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
			Object obj1 = potential2[0].createSpeciesFactory();
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
			Object obj2 = potential2[1].createSpeciesFactory();
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
	private ArrayList<String[]> removeDulplicatesFromPairSites(ArrayList<String[]> PairSitesDWithDuplicates, ParameterMapping[] Potential) {
		
		
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
	
	
	
	private ArrayList<String[]> removeSimilarCrossPairSites(ArrayList<String[]> PairSitesDWithoutDuplicates, ParameterMapping[] Potential) {
		
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
