/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;

import etomica.api.IPotential;
import etomica.api.IPotentialMolecular;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.potential.P2LJQ;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.virial.GUI.models.EnumSiteName;
import etomica.virial.GUI.models.IMolecularModel_SpeciesFactory;
import etomica.virial.GUI.models.ModelSelectedSpecies;
import etomica.virial.GUI.models.ModelSimulationEnvironment;
import etomica.virial.GUI.views.ViewAlertMsgBox;
import etomica.virial.MCMoveClusterTorsionMulti;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

public class BuilderCollectionPotential {

	private ViewAlertMsgBox messageAlert;
	
	//Main Flag to continue Simulation after potentials match
	private boolean potentialCheckFlag = false;
	
	private boolean isSpecies1AtomicFlag;
	private boolean isSpecies2AtomicFlag;
	
	private int[] electrostaticTypeIndex;
	private char[] potentialTypeCharIndex;
	
	private boolean species1HasChargeFlag = false;
	private boolean species2HasChargeFlag = false;
	
	private boolean species1HasMomentFlag = false;
	private boolean species2HasMomentFlag = false;
	
	
	
	private int[][] speciesInteractionIndex;

	
	private ModelSimulationEnvironment simEnvironment;
	private ModelSelectedSpecies speciesDataModel;
	private IMolecularModel_SpeciesFactory[] speciesDM;
	
	
	
	
	private ArrayList<String[]> siteParticipatingNBInteraction;
	private ArrayList<String[]> sitePairsWithoutDuplicates;
	private ArrayList<String[]> sitePairsWithoutDuplicateCrossPairs;
	
	//private ICollectionPotential pCollectionObject;
	private ArrayList<ICollectionPotential> pCollectionHashMap;
	
	
	@SuppressWarnings("unchecked")
	public ArrayList<ICollectionPotential> checkIfCompatible(ModelSelectedSpecies speciesDataModel, ModelSimulationEnvironment SimENV){
		
		//PObject = new PotentialObject();
		pCollectionHashMap = new ArrayList<ICollectionPotential>();
		simEnvironment = SimENV;
		this.speciesDataModel = speciesDataModel;
		int numberOfSpecies = this.speciesDataModel.getNthSpeciesAdded();
		speciesInteractionIndex = new int[numberOfSpecies*numberOfSpecies][2];
		siteParticipatingNBInteraction = new ArrayList<String[]>();
		
		
		for(int i=0;i<numberOfSpecies;i++){
			addSpeciesIndextoSites(speciesDataModel.getSpeciesDataModel(i).getPotentialSites(),i+1);
			siteParticipatingNBInteraction.add(getSitesNotParticipating(speciesDataModel.getSpeciesDataModel(i)));
		}
		
		for(int m=0;m<numberOfSpecies;m++){
			for(int n=m+1;n<numberOfSpecies;n++){
				checkForSimilarSites(siteParticipatingNBInteraction.get(m),
										siteParticipatingNBInteraction.get(n), 
											new IMolecularModel_SpeciesFactory[]{speciesDataModel.getSpeciesDataModel(m),
																					speciesDataModel.getSpeciesDataModel(n)});
			}
		}
		
		int indexCount = 0;
		for(int i=0;i<numberOfSpecies;i++){
			for(int j=0;j<numberOfSpecies;j++){
				speciesInteractionIndex[indexCount][0] = i;
				speciesInteractionIndex[indexCount][1] = j;
				if(i==j){
					if (etomica.api.IPotentialMolecular.class.isAssignableFrom(speciesDataModel.getSpeciesDataModel(i).getPotential())){
						makeMolecularPairPotentialsLike(speciesDataModel.getSpeciesDataModel(i),i);
					}else if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom(speciesDataModel.getSpeciesDataModel(i).getPotential()) && speciesDataModel.getSpeciesDataModel(i).getPotentialSites().length > 1){
						makeAtomicPairPotentialsLike(speciesDataModel.getSpeciesDataModel(i),i);
					}else if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom(speciesDataModel.getSpeciesDataModel(i).getPotential()) && speciesDataModel.getSpeciesDataModel(i).getPotentialSites().length == 1){
						makeAtomicPairPotentialsLike(speciesDataModel.getSpeciesDataModel(i),i);
					}else if(PotentialSoft.class.isAssignableFrom(speciesDataModel.getSpeciesDataModel(i).getPotential())){
						makeAtomicPairPotentialsLike(speciesDataModel.getSpeciesDataModel(i),i);
					}
				}
				indexCount++;
			}
		}
		
		for(int m=0;m<speciesInteractionIndex.length;m++){
			
			int speciesIndex1 = speciesInteractionIndex[m][0];
			int speciesIndex2 = speciesInteractionIndex[m][1];
			
			if(speciesIndex1 != speciesIndex2 && speciesIndex1 < speciesIndex2){
				makePairPotentialsUnlike(speciesIndex1,speciesIndex2);
			}
			
		}

		speciesDataModel.setSpeciesInteractionIndex(speciesInteractionIndex);
		/*
		Map<String[],ICollectionPotential> nonBondedPotentialsMap = pCollectionHashMap;
		Set<Entry<String[],ICollectionPotential>> nonBondedPotentialEntries = nonBondedPotentialsMap.entrySet();
		Iterator<Entry<String[],ICollectionPotential>> nonBondedPotentialsItr = nonBondedPotentialEntries.iterator();
		
		while(nonBondedPotentialsItr.hasNext()){
			Map.Entry nonBondedPotentialsEntry = (Map.Entry) nonBondedPotentialsItr.next();
			String[] nonBondedPotentialsMapKey= (String[]) nonBondedPotentialsEntry.getKey();
			System.out.println(nonBondedPotentialsMapKey[0]+" "+nonBondedPotentialsMapKey[1]+" "+ nonBondedPotentialsEntry.toString()+"\n");
		}*/
		
		if(simEnvironment.isStopSimulation()){
			return null;
		}else{
			return pCollectionHashMap;
		}
		
		//System.out.println("All like Pairs created :D");
		
		
		
		//return null;
	}
	
	
	public void makePairPotentialsUnlike(int speciesIndex1, int speciesIndex2){
		speciesDM = new IMolecularModel_SpeciesFactory[2];
		speciesDM[0] = speciesDataModel.getSpeciesDataModel(speciesIndex1);
		speciesDM[1] = speciesDataModel.getSpeciesDataModel(speciesIndex2);
		
		potentialTypeCharIndex = new char[2];
		electrostaticTypeIndex = new int[2];
		
		for(int i = 0;i < 2;i++){
			//We figure details abt each of the potential
			
			//If molecular or atomic
			if (etomica.api.IPotentialMolecular.class.isAssignableFrom(speciesDM[i].getPotential())){

				//If potentialsites existing are greater than 1, although we dont have a intrabonded or intra non-bonded potential, but we have 
				//cross potentials
				if(i==0){
					setSpecies1AtomicFlag(false);
				}
				if(i==1){
					setSpecies2AtomicFlag(false);
					
				}
				potentialTypeCharIndex[i] = 'M';
				//PObject.setpIntraGroupII(null, i);
				//pIntraTargetGroup[i] = null;
			}
			if(etomica.potential.Potential2SoftSpherical.class.isAssignableFrom(speciesDM[i].getPotential()) || etomica.potential.PotentialSoft.class.isAssignableFrom(speciesDM[i].getPotential())){
				
					//Inter-bonded potentials is to be gathered especially for alkanes, alcohol, etc
					if(i==0){
						setSpecies1AtomicFlag(true);
						potentialTypeCharIndex[0] = 'A';
						
					}
					if(i==1){
						setSpecies2AtomicFlag(true);
						potentialTypeCharIndex[1] = 'A';
						
					}
					
			}
			
			//If potential describing the interaction incluses charges or moments
			String[] tempArray = speciesDM[i].getParametersArray();
			for(int j=0; j< tempArray.length;j++){
				if(tempArray[j].equals("CHARGE")){
					if(i == 0){
						species1HasChargeFlag = true;}
					if(i == 1){
						species2HasChargeFlag = true;}
					
				}
				if(tempArray[j].equals("MOMENT")||tempArray[j].equals("MOMENTSQR")){
					//Dipole /Quadrapole moment exists!
					if(i == 0){
						species1HasMomentFlag = true;}
					if(i == 1){
						species2HasMomentFlag = true;}

				}
			}
		}
		if(potentialTypeCharIndex[0] == potentialTypeCharIndex[1]){
			//Both are molecular Potentials or Both are Atomic Potentials
			if(potentialTypeCharIndex[0] == 'M'){
				if(speciesDM[0].getPotential().equals(speciesDM[1].getPotential())){
					potentialCheckFlag = true;	
					makeMolecularPairPotentialsUnlike(speciesDM[0], speciesDM[1],speciesIndex1,speciesIndex2);
					
				}
				else if( speciesDM[0].getNonBondedInteractionModel().equals(speciesDM[1].getNonBondedInteractionModel()) ){
					//CrossPotentials are site-site interactions and modelled using either lennard Jones or Exp-6 potential
					potentialCheckFlag = true;
					
					makeAtomicPairPotentialsUnlike(speciesDM[0], speciesDM[1],speciesIndex1,speciesIndex2,false);
				}
			}
			
			
			if(potentialTypeCharIndex[0] == 'A'){
				
				if(speciesDM[0].getPotential().equals(speciesDM[1].getPotential())){
					potentialCheckFlag = true;
					makeAtomicPairPotentialsUnlike(speciesDM[0], speciesDM[1],speciesIndex1,speciesIndex2,true);
				

				}
				else if( speciesDM[0].getNonBondedInteractionModel().equals(speciesDM[1].getNonBondedInteractionModel()) ){
					potentialCheckFlag = true;
					makeAtomicPairPotentialsUnlike(speciesDM[0], speciesDM[1],speciesIndex1,speciesIndex2,false);
					
				}
			}
		}
		else{
			if( speciesDM[0].getNonBondedInteractionModel().equals(speciesDM[1].getNonBondedInteractionModel()) ){
				potentialCheckFlag = true;
				makeAtomicPairPotentialsUnlike(speciesDM[0], speciesDM[1],speciesIndex1,speciesIndex2,false);
			}else{
				potentialCheckFlag = false;
			}
			
		}
		
		if(potentialCheckFlag){
			//Condition for electrostatic interaction included
			electrostaticTypeIndex = getElectrostatic();
			if(electrostaticTypeIndex[0] == electrostaticTypeIndex[1] || electrostaticTypeIndex[0] == 0 || electrostaticTypeIndex[1] == 0){
				System.out.println("Run");
				
			}
			else{
				//stop simulation
				System.out.println("Stop");
				simEnvironment.setStopSimulation(true);
				
			}
		}
		else{
			//stop simulation
			System.out.println("Stop");
			simEnvironment.setStopSimulation(true);
		}
		
	}
	
	
	
	private void makeAtomicPairPotentialsUnlike(IMolecularModel_SpeciesFactory species1DataModel, IMolecularModel_SpeciesFactory species2DataModel,
													int species1Index,int species2Index,boolean samePotentialFlag){

		ArrayList<String[]> TempPairSites = makePairSites(siteParticipatingNBInteraction.get(species1Index),siteParticipatingNBInteraction.get(species2Index));
		
		sitePairsWithoutDuplicates = removeDulplicatesFromPairSites(TempPairSites);
	
		//Removes Similar Cross Pairs Example : C,H and H,C
		sitePairsWithoutDuplicateCrossPairs = removeSimilarCrossPairSites(sitePairsWithoutDuplicates);
	
		IMolecularModel_SpeciesFactory[] PotentialTemp = new IMolecularModel_SpeciesFactory[2];
		PotentialTemp[0] = species1DataModel;
		PotentialTemp[1] = species2DataModel;
		
		
		ISpecies[] tempSpecies = new ISpecies[2];
		tempSpecies[0] = speciesDataModel.getSpecies(species1Index);
		tempSpecies[1] = speciesDataModel.getSpecies(species2Index);
		
		
		CollectionPotentialAtomicUnlike collectionPotentialAtomicUnlike = new CollectionPotentialAtomicUnlike(species1Index,species2Index);
		collectionPotentialAtomicUnlike.setHashMapAtomTypesUnlikePairs(getPairAtomTypes(sitePairsWithoutDuplicates,tempSpecies));
		
		if(samePotentialFlag == true){
			if(species1DataModel.getPotential().equals(P2LennardJones.class)){
				//PotentialSets = getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class);
				getPairSitePotentials(sitePairsWithoutDuplicateCrossPairs,PotentialTemp, P2LennardJones.class,
						collectionPotentialAtomicUnlike);
			}
			if(species1DataModel.getPotential().equals(P2LJQ.class)){
				//PotentialSets = getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class);
				getPairSitePotentials(sitePairsWithoutDuplicateCrossPairs,PotentialTemp, P2LJQ.class,
						collectionPotentialAtomicUnlike);
			}
		}else{
			if(species1DataModel.getNonBondedInteractionModel() == "LennardJones"){
				getPairSitePotentials(sitePairsWithoutDuplicateCrossPairs,PotentialTemp, P2LennardJones.class,
						collectionPotentialAtomicUnlike);
			}
			
			if(species1DataModel.getNonBondedInteractionModel() == "LennardJonesWithQuadrapole"){
				//PotentialSets = getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class);
				getPairSitePotentials(sitePairsWithoutDuplicateCrossPairs,PotentialTemp, P2LJQ.class,
						collectionPotentialAtomicUnlike);
			}
		}
		
		
		pCollectionHashMap.add(collectionPotentialAtomicUnlike);
		
	}
	
	private void makeAtomicPairPotentialsLike(IMolecularModel_SpeciesFactory species1DataModel, int speciesIndex){

		// Make Pair Sites
		ArrayList<String[]> TempPairSites = makePairSites(siteParticipatingNBInteraction.get(speciesIndex),siteParticipatingNBInteraction.get(speciesIndex));
	
		sitePairsWithoutDuplicates = removeDulplicatesFromPairSites(TempPairSites);
	
		//Removes Similar Cross Pairs Example : C,H and H,C
		sitePairsWithoutDuplicateCrossPairs = removeSimilarCrossPairSites(sitePairsWithoutDuplicates);
	
		IMolecularModel_SpeciesFactory[] PotentialTemp = new IMolecularModel_SpeciesFactory[2];
		PotentialTemp[0] = species1DataModel;
		
		ISpecies[] tempSpecies = new ISpecies[2];
		tempSpecies[0] = this.speciesDataModel.getSpecies(speciesIndex);
		
		CollectionPotentialAtomicLike collectionPotentialAtomicLike = new CollectionPotentialAtomicLike(speciesIndex);
		collectionPotentialAtomicLike.setHashMapAtomTypeLikePairs(getPairAtomTypes(sitePairsWithoutDuplicates,tempSpecies));
		if(species1DataModel.getNonBondedInteractionModel() == "LennardJones"){
			//PotentialSets = getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class);
			getPairSitePotentials(sitePairsWithoutDuplicateCrossPairs,PotentialTemp, P2LennardJones.class,
									collectionPotentialAtomicLike);
		}
		if(species1DataModel.getNonBondedInteractionModel() == "LennardJonesWithQuadrapole"){
			//PotentialSets = getPairSitePotentials(UnlikePairsWithoutCrossPairs,potential, P2LennardJones.class);
			getPairSitePotentials(sitePairsWithoutDuplicateCrossPairs,PotentialTemp, P2LJQ.class,
									collectionPotentialAtomicLike);
		}
		
		
		if(simEnvironment.getLengthAlkaneChain(speciesIndex) > 2){
			ACollectionBondedPotential bInteractions = new CollectionBondedPotentialAlkanes(simEnvironment.getLengthAlkaneChain(speciesIndex));
			bInteractions.addBondedPotentialSets(collectionPotentialAtomicLike,species1DataModel.getSpace(),1);
			if(simEnvironment.getLengthAlkaneChain(speciesIndex) > 3){
				collectionPotentialAtomicLike.setMCMoveClusterTorsionMulti(new MCMoveClusterTorsionMulti[2]);
			}
		}
		else if(species1DataModel.getClass().getName().contains("Methanol") || species1DataModel.getClass().getName().contains("Ethanol")){
			ACollectionBondedPotential bInteractions = new CollectionBondedPotentialAlcohols();
			bInteractions.addBondedPotentialSets(collectionPotentialAtomicLike,species1DataModel.getSpace(),1);
		}
		
		pCollectionHashMap.add(collectionPotentialAtomicLike);
		
		
	}
	
	private void makeMolecularPairPotentialsLike(IMolecularModel_SpeciesFactory species1DataModel, int speciesIndex){
		boolean OneConstructor = hasSingleConstructorForPotential(species1DataModel.getPotential());
		int numberofParameters = species1DataModel.getParametersArray().length;
		
	
		@SuppressWarnings("rawtypes")
		Class[] ParamClass = new Class[numberofParameters+1];
		ParamClass[0] = Space.class;
		for(int j = 1;j<=numberofParameters;j++){
			ParamClass[j] = Double.TYPE;
		}
	
		Object[] ParamValueObj = new Object[numberofParameters+1];
		ParamValueObj[0] = species1DataModel.getSpace();
		if(!OneConstructor){
			for(int j=0;j<species1DataModel.getParametersArray().length;j++){
				ParamValueObj[j+1]=species1DataModel.getDoubleDefaultParameters(species1DataModel.getParametersArray()[j].toUpperCase()+
																					species1DataModel.getPotentialSiteAtIndex(0));
			}
		}
		CollectionPotentialMolecularLike collectionPotentialMolecularLike = new CollectionPotentialMolecularLike(speciesIndex);
		try{
				if(OneConstructor){
					//Species1Molecular = potential[0].getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]);
					collectionPotentialMolecularLike.setPotentialMolecularNonBondedLike((IPotentialMolecular) species1DataModel.getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]));
					
				}else{
					collectionPotentialMolecularLike.setPotentialMolecularNonBondedLike((IPotentialMolecular) species1DataModel.getPotential().getConstructor(ParamClass).newInstance(ParamValueObj));
				}
		}
		catch(Exception E){
			E.printStackTrace();
		}
		
		pCollectionHashMap.add(collectionPotentialMolecularLike);
		
	
	}
	
	private void makeMolecularPairPotentialsUnlike(IMolecularModel_SpeciesFactory species1DataModel, IMolecularModel_SpeciesFactory species2DataModel,
			int species1Index,int species2Index){
		
		boolean OneConstructor = hasSingleConstructorForPotential(species1DataModel.getPotential());
		int numberofParameters = species1DataModel.getParametersArray().length;
		Object[][] ParamValueCrossObj = new Object[numberofParameters][3];
		Object[] ParamValueObj = new Object[numberofParameters+1];
		@SuppressWarnings("rawtypes")
		Class[] ParamClass = new Class[numberofParameters+1];
		ParamClass[0] = Space.class;
		for(int j = 1;j<=numberofParameters;j++){
			ParamClass[j] = Double.TYPE;
		}
		
		ParamValueObj[0] = species1DataModel.getSpace();
		
		if(!OneConstructor && species1DataModel.getPotentialSites().length == 1 ){
				for(int j=0;j<species1DataModel.getParametersArray().length;j++){
					ParamValueCrossObj[0][j]=species1DataModel.getParametersArray()[j].toUpperCase();
					ParamValueCrossObj[1][j]=species1DataModel.getDoubleDefaultParameters(species1DataModel.getParametersArray()[j].toUpperCase()+species1DataModel.getPotentialSiteAtIndex(0));
				}
		}
		
		if(!OneConstructor && species2DataModel.getPotentialSites().length == 1 ){
			for(int j=0;j<species2DataModel.getParametersArray().length;j++){
				ParamValueCrossObj[2][j]=species2DataModel.getDoubleDefaultParameters(species2DataModel.getParametersArray()[j].toUpperCase()+species2DataModel.getPotentialSiteAtIndex(0));
			}
		}
		
		String valueA = null;
		String valueB = null;
		if(!OneConstructor){
			for(int i=0;i < ParamValueCrossObj.length;i++){
				valueA = ParamValueCrossObj[1][i].toString();
				double ValueA = Double.parseDouble(valueA);
				valueB = ParamValueCrossObj[2][i].toString();
				double ValueB = Double.parseDouble(valueB);		
				ParamValueObj[i+1]=GetMixingRules(ValueA,ValueB,ParamValueCrossObj[0][i].toString());
			}
		}
		
		CollectionPotentialMolecularUnlike collectionPotentialMolecularUnlike = new CollectionPotentialMolecularUnlike(species1Index,species2Index);
		try {
			if(OneConstructor){
				collectionPotentialMolecularUnlike.setPotentialMolecularUnlike((IPotentialMolecular) species1DataModel.getPotential().getConstructor(ParamClass[0]).newInstance(ParamValueObj[0]));
			}else{
				collectionPotentialMolecularUnlike.setPotentialMolecularUnlike((IPotentialMolecular) species1DataModel.getPotential().getConstructor(ParamClass).newInstance(ParamValueObj));
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		pCollectionHashMap.add(collectionPotentialMolecularUnlike);
		
	}
	
	public void addSpeciesIndextoSites(String[] sitesOnSpecies, int speciesIndex){
		for(int i=0;i<sitesOnSpecies.length;i++){
			sitesOnSpecies[i] = sitesOnSpecies[i] + "-" + Integer.toString(speciesIndex);
		}
	}
	
	
	public void checkForSimilarSites(String[] sitesOnSpecies1,String[] sitesOnSpecies2, IMolecularModel_SpeciesFactory[] speciesSelected){
		
		
		
		for(int i=0;i<sitesOnSpecies1.length;i++){
			for(int j=0;j<sitesOnSpecies2.length;j++){
				StringTokenizer dataString1 = new StringTokenizer(sitesOnSpecies1[i],"-");
				StringTokenizer dataString2 = new StringTokenizer(sitesOnSpecies2[j],"-");
				String potentialSite1 = dataString1.nextToken().trim();
				String potentialSite2 = dataString2.nextToken().trim();
				if(potentialSite1.equals(potentialSite2)){
					if(getPairSiteValidation(potentialSite2,speciesSelected)){
						//newP1[i] = newP1[i] + "-a";
						//newP2[j] = newP2[j] + "-b";
						for(int k=0;k<speciesSelected[1].getPotentialSites().length;k++){
							if(speciesSelected[1].getPotentialSiteAtIndex(k)==sitesOnSpecies2[j]){
								speciesSelected[1].getPotentialSites()[k] = sitesOnSpecies1[i];
							}
						}
						sitesOnSpecies2[j] = sitesOnSpecies1[i];
						
						 
					}
				}
			}
		}
	}
	
	
	/*
	 * Function Retrieves the IAtomTypes for Pair Sites
	 */
	
	
	private void getSameSitesOnBothSpecies(String[] newP1, String[] newP2,
			IMolecularModel_SpeciesFactory species1DataModel,IMolecularModel_SpeciesFactory species2DataModel) {
		IMolecularModel_SpeciesFactory[] speciesDataModel = new IMolecularModel_SpeciesFactory[2];
		speciesDataModel[0] = species1DataModel;
		speciesDataModel[1] = species2DataModel;
		
		for(int i=0;i<newP1.length;i++){
			for(int j=0;j<newP2.length;j++){
				if(newP1[i]==newP2[j]){
					if(!getPairSiteValidation(newP1[i],speciesDataModel)){
						//newP1[i] = newP1[i] + "-a";
						//newP2[j] = newP2[j] + "-b";
					}else{
						newP2[j] = newP1[i];
						 
					}
				}
			}
		}
		
		
		
	}



	@SuppressWarnings({ "unchecked", "rawtypes" })
	private void getPairSitePotentials(
			ArrayList<String[]> unlikePairsWithoutCrossPairs,
			IMolecularModel_SpeciesFactory[] potential,Class TargetPotential,ICollectionPotential potentialCollection) {
		String[] ParametersForPotential = null;
		if(TargetPotential.equals(P2LennardJones.class)){
			ParametersForPotential = new String[]{"SIGMA","EPSILON"};
			hasSingleConstructorForPotential(TargetPotential);
		}
		if(TargetPotential.equals(P2LJQ.class)){
			ParametersForPotential = new String[]{"SIGMA","EPSILON","MOMENTSQR"};
			hasSingleConstructorForPotential(TargetPotential);
		}
		Double[] Site12 = new Double[ParametersForPotential.length];
		
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
			ParamClass[0] = Space.class;
			for(int j1 = 1;j1<=ParametersForPotential.length;j1++){
				ParamClass[j1] = Double.TYPE;
			}
			
			Object[] ParamValueObj = new Object[ParametersForPotential.length+1];
			ParamValueObj[0] = potential[0].getSpace();
			for(int m=0;m<Site12.length;m++){
				ParamValueObj[m+1]=Site12[m];
			}
			
			try {
				//PairPotentialSets.put(unlikePairsWithoutCrossPairs.get(i),TargetPotential.getConstructor(ParamClass).newInstance(ParamValueObj));
				//PairPotentialSets.addToHashMapPotentialNonBonded(unlikePairsWithoutCrossPairs.get(i),TargetPotential.getConstructor(ParamClass).newInstance(ParamValueObj));
				if(potentialCollection instanceof CollectionPotentialAtomicLike){
					((CollectionPotentialAtomicLike) potentialCollection).
							addToHashMapPotentialNonBonded(unlikePairsWithoutCrossPairs.get(i), (IPotential) TargetPotential.getConstructor(ParamClass).newInstance(ParamValueObj));
				}else if(potentialCollection instanceof CollectionPotentialAtomicUnlike){
					((CollectionPotentialAtomicUnlike) potentialCollection).
					addToHashMapPotentialNonBonded(unlikePairsWithoutCrossPairs.get(i), (IPotential) TargetPotential.getConstructor(ParamClass).newInstance(ParamValueObj));
				}
			} catch (Exception e) {
				
				e.printStackTrace();
			}
		}
		
		
		
		//return PairPotentialSets;
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
    private HashMap<String[], AtomType[]> getPairAtomTypes(ArrayList<String[]> unLikePairs,
                                                           ISpecies[] species) {

		ISpecies species1 = species[0];		
		ISpecies species2;
		
		if(species[1] != null){
			species2 = species[1];
		}
		else{
			species2 = species1;}
		String Site1Name = "";
		String Site2Name = "";
        AtomType Site1Atom;
        AtomType Site2Atom;
        HashMap<String[], AtomType[]> AtomSets = new HashMap<String[], AtomType[]>();

		for(int i=0;i<unLikePairs.size();i++){
			String Site1 = unLikePairs.get(i)[0];
			String Site2 = unLikePairs.get(i)[1];
			//if(Site1.contains("-a") || Site1.contains("-b")){
			if(Site1.contains("-")){
				String[] SiteA = Site1.split("-");
				Site1 = SiteA[0];
			}
			//if(Site2.contains("-a") || Site2.contains("-b")){
			if(Site2.contains("-")){
				String[] SiteB = Site2.split("-");
				Site2 = SiteB[0];
			}
            AtomType[] SiteAtoms = new AtomType[2];
            for(EnumSiteName site: EnumSiteName.values()){
				if(site.toString().toUpperCase().equals(Site1.toUpperCase())){
					Site1Name = site.getSite();
				}
				if(site.toString().toUpperCase().equals(Site2.toUpperCase())){
					Site2Name = site.getSite();
				}
			}

			Method[] MethodS1 = species1.getClass().getMethods();
			
			
			for(int j=0;j<MethodS1.length;j++){
                if (MethodS1[j].getReturnType().equals(AtomType.class)) {
                    if(MethodS1[j].toString().toUpperCase().contains(Site1.toUpperCase()+"TYPE")|| MethodS1[j].toString().toUpperCase().contains(Site1Name.toUpperCase()+"TYPE")){
						try {
                            Site1Atom = (AtomType) MethodS1[j].invoke(species1);
                            SiteAtoms[0] = Site1Atom;
						} catch (Exception e) {
							e.printStackTrace();
						}
						
					}
				}
			}
			
			Method[] MethodS2 = species2.getClass().getMethods();
			
			for(int k=0;k<MethodS2.length;k++){
                if (MethodS2[k].getReturnType().equals(AtomType.class)) {
                    if(MethodS2[k].toString().toUpperCase().contains(Site2.toUpperCase()+"TYPE")|| MethodS2[k].toString().toUpperCase().contains(Site2Name.toUpperCase()+"TYPE")){
						System.out.println(Site2+" "+ MethodS2[k].toString());
						try {
                            Site2Atom = (AtomType) MethodS2[k].invoke(species2);
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
		//String[][] PairSitesDraft2on1 = new String[newP2.length*newP1.length][2];
		
		
		int index= 0;
		
		for(int i=0;i<newP1.length;i++){
			for(int j=0;j<newP2.length;j++){
				PairSitesDraft1on2[index][0]=newP1[i];
				PairSitesDraft1on2[index][1]=newP2[j];
				index++;
			}
		}
		/*for(int i=0;i<newP2.length;i++){
			for(int j=0;j<newP1.length;j++){
				PairSitesDraft2on1[index1][0]=newP2[i];
				PairSitesDraft2on1[index1][1]=newP1[j];
				index1++;
			}
		}*/
		
		ArrayList<String[]> PairSitesD = new ArrayList<String[]>();
		for(int i=0;i<PairSitesDraft1on2.length;i++){
			PairSitesD.add(PairSitesDraft1on2[i]);
		}
		/*for(int i=PairSitesDraft1on2.length;i<PairSitesDraft1on2.length+PairSitesDraft2on1.length;i++){
			PairSitesD.add(PairSitesDraft2on1[i - PairSitesDraft1on2.length]);
		}*/
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
		if((species1HasChargeFlag && species2HasChargeFlag)){
			//run simulation
			
			return new int[]{1,1};
		}
		else if(species1HasMomentFlag && species2HasMomentFlag){
			//makeMolecularPairPotentials(potential);
			return new int[]{2,2};
		}
		else if(!species1HasMomentFlag && !species2HasMomentFlag && !species1HasChargeFlag && !species2HasChargeFlag){
			return new int[]{0,0};
		}
		//-----------------No charge or moment on species 1 but charge or moment on species 2
		
		else if(!species1HasChargeFlag && !species1HasMomentFlag && species2HasMomentFlag){
			return new int[]{0,2};
		}
		
		else if(!species1HasChargeFlag && !species1HasMomentFlag && species2HasChargeFlag ){
			return new int[]{0,1};
		}
		
		//------------------Charge or moment on species 1 but no charge or moment on species 2
		else if(species1HasMomentFlag && !species2HasChargeFlag && !species2HasMomentFlag){
			return new int[]{2,0};
		}
		
		else if(species1HasChargeFlag & !species2HasChargeFlag && !species2HasMomentFlag){
			return new int[]{1,0};
		}
		//--------------------------------Charge on species 1 and Moment on species 2 || Moment on species 1 and charge on species 2
		else if(species1HasChargeFlag && species2HasMomentFlag){
			return new int[]{1,2};
		}
		else if(species1HasMomentFlag && species2HasChargeFlag){
			return new int[]{2,1};
		}
		else{
		return null;
		}
	}
	
	
	public String[] getSitesNotParticipating(IMolecularModel_SpeciesFactory Potential){
		ArrayList<String> TempArray = new ArrayList<String>();
		String[] NewP;
		for(int i=0;i < Potential.getPotentialSites().length;i++){
			//StringTokenizer potentialSiteString = null;
			
			String PotentialSite = Potential.getPotentialSites()[i];
			//potentialSiteString = new StringTokenizer((String)PotentialSite,"-");
			for(int j=0;j<Potential.getParametersArray().length;j++){
				String Param = Potential.getParametersArray()[j];
				Double Value = Potential.getDoubleDefaultParameters(Param+PotentialSite);
				//Double Value = Potential.getDoubleDefaultParameters(Param+PotentialSite);
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
	
	
	
	
	
	
	public boolean getPairSiteValidation(String PotentialSite, IMolecularModel_SpeciesFactory[] Potential){
		
			Object[][] ParamValueCrossObj = null;
			int ParamCheck = 0;
		
			if(Potential[0].getNonBondedInteractionModel() =="LennardJones"){
				ParamValueCrossObj = new Object[2][3];
				ParamCheck=2;
			}
			if(Potential[0].getNonBondedInteractionModel() =="LennardJonesWithQuadrapole"){
				ParamValueCrossObj = new Object[3][3];
				ParamCheck=3;
			}
			
			for(int k=0;k<2;k++){
				for(int l=0;l<Potential[k].getPotentialSites().length;l++){
					StringTokenizer dataString1 = new StringTokenizer(Potential[k].getPotentialSiteAtIndex(l),"-");
					String potentialSite = dataString1.nextToken().trim();
					if(PotentialSite.equals(potentialSite)){
						for(int m=0;m<Potential[k].getParametersArray().length;m++){
							String Param = Potential[k].getParametersArray()[m];
							if(Param == "SIGMA"|| Param == "EPSILON"|| Param == "MOMENTSQR" ){
								ParamValueCrossObj[m][0] = Param;
								if(k==0){
									ParamValueCrossObj[m][1] = Potential[k].getDoubleDefaultParameters(Param+Potential[k].getPotentialSiteAtIndex(l));}
								if(k==1){
									ParamValueCrossObj[m][2] = Potential[k].getDoubleDefaultParameters(Param+Potential[k].getPotentialSiteAtIndex(l));}
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

        return ConfirmIndexFlag == ParamCheck;
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
		return isSpecies1AtomicFlag;
	}



	public void setSpecies1AtomicFlag(boolean species1AtomicFlag) {
		isSpecies1AtomicFlag = species1AtomicFlag;
	}



	public boolean isSpecies2AtomicFlag() {
		return isSpecies2AtomicFlag;
	}



	public void setSpecies2AtomicFlag(boolean species2AtomicFlag) {
		isSpecies2AtomicFlag = species2AtomicFlag;
	}


	public boolean isPotentialSite(String PotentialSite,IMolecularModel_SpeciesFactory Potential){
		
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
