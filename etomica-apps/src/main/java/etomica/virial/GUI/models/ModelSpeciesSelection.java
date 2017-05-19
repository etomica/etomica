/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;

import javax.swing.ListSelectionModel;

public class ModelSpeciesSelection {
	
	private static final String[] intialDisplayList = {"---No species selected---","Press \'Add\' to select a species"};
	private static final String[] intialPotentialDisplayList = {"---No potentials selected---","Press \'Add\' to select a potential"};
	
	private static final String[] species_List = {"LJ","CO2","Methane","Ethane","Propane","Higher n-Alkanes","Naphthalene","Water"};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesLJ_PotentialsList = {MolecularModelLJ_SpeciesLJ.class,MolecularModelLJQ_SpeciesLJ.class,MolecularModel2CLJQ_SpeciesLJ.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesCO2_PotentialsList = {MolecularModel2CLJQ_SpeciesCO2.class,MolecularModelEMP2_SpeciesCO2.class,MolecularModelTrappe_SpeciesCO2.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesMethane_PotentialList = {MolecularModelTRAPPE_SpeciesAlkane.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesEthane_PotentialList = {MolecularModelTRAPPE_SpeciesAlkane.class, MolecularModelSKS_SpeciesAlkane.class, MolecularModel2CLJQ_SpeciesEthane.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesPropane_PotentialList = {MolecularModelTRAPPE_SpeciesAlkane.class,MolecularModelSKS_SpeciesAlkane.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesNAlkane_PotentialList = {MolecularModelTRAPPE_SpeciesAlkane.class, MolecularModelSKS_SpeciesAlkane.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesH2O_PotentialList = {MolecularModelSPCE_SpeciesH2O.class};
	
	private String[] speciesIdentifier;
	
	//AddFlag in Species1
	private boolean[] speciesAlreadyAdded;

	
	//ChangesSpecies in Species1
	private boolean speciesSelectionChanged;
	
	
	//Indexs for storing JList.selectedIndex()
	private int nthSpeciesSelectedIndex;
	private int[] speciesSelectedIndex;
	private int[] potentialsSelectedIndex;

	//Species-sites
	private String[] sitesOnSpecies;

	private ArrayList<String> moleculeListDisplayString;
	
	private String[] potentialsStringList;
	private Class[] potentialsClassList;
	
	private String[] countOfMoleculesStrings = {"0","1","2","3","4","5","6","7","8"};

	
	private boolean removeSpeciesAfterAlertFlag;
	//potential Parameters;
	private String[] potentialParamString;

	private String[][] potentialParamAndValuesString;

	private IMolecularModel_SpeciesFactory selectedSpecies;
	
	private Object[] constructorObject;

	public ModelSpeciesSelection(){
		reset();
	}
	
	public void reset() {
		speciesIdentifier = new String[10];
		speciesSelectedIndex = new int[10];
		potentialsSelectedIndex = new int[10];
		speciesAlreadyAdded = new boolean[10];

		for (int i=0;i <10;i++){
			speciesAlreadyAdded[i]=false;
		}
		speciesSelectionChanged = false;
		removeSpeciesAfterAlertFlag = false;
		selectedSpecies = null;
		
	}
	
	public Object[] getConstructorObject() {
		return constructorObject;
	}

	public void setConstructorObject(Integer constructorObject) {
		this.constructorObject = new Object[1];
		this.constructorObject[0] = constructorObject;
	}

	public String getSpeciesIdentifier(int index) {
		return speciesIdentifier[index];
	}

	public void setSpeciesIdentifier(String speciesIdentifier, int index) {
		this.speciesIdentifier[index] = speciesIdentifier;
	}

	public String[] getPotentialsList() {
		return potentialsStringList;
	}
	
	
	
	public int getNthSpeciesSelectedIndex() {
		return nthSpeciesSelectedIndex;
	}

	public void setNthSpeciesSelectedIndex(int nthSpeciesSelectedIndex) {
		this.nthSpeciesSelectedIndex = nthSpeciesSelectedIndex;
	}

	public Class[] getPotentialsClassList() {
		return potentialsClassList;
	}

	public void setPotentialsClassList(Class[] potentialsClassList) {
		this.potentialsClassList = potentialsClassList;
	}

	public int[] getSpeciesSelectedIndex() {
		return speciesSelectedIndex;
	}
	
	public int getSpeciesSelectedIndex(int index) {
		return speciesSelectedIndex[index];
	}

	public void setSpeciesSelectedIndex(int[] speciesSelectedIndex) {
		this.speciesSelectedIndex = speciesSelectedIndex;
	}
	
	public void setSpeciesSelectedIndex(int speciesSelectedIndex, int index) {
		this.speciesSelectedIndex[index] = speciesSelectedIndex;
	}

	public int[] getPotentialsSelectedIndex() {
		return potentialsSelectedIndex;
	}

	public void setPotentialsSelectedIndex(int[] potentialsSelectedIndex) {
		this.potentialsSelectedIndex = potentialsSelectedIndex;
	}
	
	public void setPotentialsSelectedIndex(int potentialsSelectedIndex, int index) {
		this.potentialsSelectedIndex[index] = potentialsSelectedIndex;
	}

	
	public static String[] getIntialdisplaylist() {
		return intialDisplayList;
	}

	public static String[] getIntialpotentialdisplaylist() {
		return intialPotentialDisplayList;
	}

	public static String[] getSpeciesList() {
		return species_List;
	}

	public static Class[] getSpeciesljPotentialslist() {
		return speciesLJ_PotentialsList;
	}

	public static Class[] getSpeciesco2Potentialslist() {
		return speciesCO2_PotentialsList;
	}

	public static Class[] getSpeciesmethanePotentiallist() {
		return speciesMethane_PotentialList;
	}

	public static Class[] getSpeciesethanePotentiallist() {
		return speciesEthane_PotentialList;
	}

	public static Class[] getSpeciespropanePotentiallist() {
		return speciesPropane_PotentialList;
	}

	public static Class[] getSpeciesnalkanePotentiallist() {
		return speciesNAlkane_PotentialList;
	}
	
	public static Class[] getSpeciesh2oPotentiallist() {
		return speciesH2O_PotentialList;
	}
	
	public ArrayList<String> getMoleculeListDisplayString() {
		return moleculeListDisplayString;
	}

	public void setMoleculeListDisplayString(
			ArrayList<String> moleculeListDisplayString) {
		this.moleculeListDisplayString = moleculeListDisplayString;
	}
	
	
	public boolean isSpeciesAlreadyAdded(int index) {
		return speciesAlreadyAdded[index];
	}

	public void setSpeciesAlreadyAdded(boolean species1AlreadyAdded, int index) {
		this.speciesAlreadyAdded[index] = species1AlreadyAdded;
	}
	
	
	public boolean isSpeciesSelectionChanged() {
		return speciesSelectionChanged;
	}

	public void setSpeciesSelectionChanged(boolean speciesSelectionChanged) {
		this.speciesSelectionChanged = speciesSelectionChanged;
	}
	
	

	public IMolecularModel_SpeciesFactory getSelectedSpecies() {
		return selectedSpecies;
	}

	public void setSelectedSpecies(IMolecularModel_SpeciesFactory selectedSpecies) {
		this.selectedSpecies = selectedSpecies;
	}

	public String[] getSitesOnSpecies() {
		return sitesOnSpecies;
	}

	public void setSitesOnSpecies(String[] sitesOnSpecies) {
		this.sitesOnSpecies = sitesOnSpecies;
	}

	public String[] getPotentialParamString() {
		return potentialParamString;
	}

	public void setPotentialParamString(String[] potentialParamString) {
		this.potentialParamString = potentialParamString;
	}

	public String[][] getPotentialParamAndValuesString() {
		return potentialParamAndValuesString;
	}

	public void setPotentialParamAndValuesString(
			String[][] potentialParamAndValuesString) {
		this.potentialParamAndValuesString = potentialParamAndValuesString;
	}

	
	public String[] getCountOfMoleculesStrings() {
		return countOfMoleculesStrings;
	}

	public void setCountOfMoleculesStrings(String[] countOfMoleculesStrings) {
		this.countOfMoleculesStrings = countOfMoleculesStrings;
	}
	

	public boolean isRemoveSpeciesAfterAlertFlag() {
		return removeSpeciesAfterAlertFlag;
	}

	public void setRemoveSpeciesAfterAlertFlag(boolean removeSpeciesAfterAlertFlag) {
		this.removeSpeciesAfterAlertFlag = removeSpeciesAfterAlertFlag;
	}

	public void setPotentialList(Class[] potentialList){
		this.potentialsClassList = potentialList;
		this.potentialsStringList = new String[potentialList.length];
        for (int i=0; i<potentialsStringList.length; i++) {
		   try {
			   try {
				   potentialsStringList[i] = (String)potentialList[i].getMethod("getMolecularModelDisplayName", new Class[0]).invoke(potentialList[i].getConstructor().newInstance());
					} 
			   catch (InstantiationException e1) {
				   // TODO Auto-generated catch block
				   e1.printStackTrace();
			   }
		   } catch (IllegalArgumentException e1) {
			   // TODO Auto-generated catch block
			   e1.printStackTrace();
		   } catch (SecurityException e1) {
			   // TODO Auto-generated catch block
			   e1.printStackTrace();
		   } catch (IllegalAccessException e1) {
			   // TODO Auto-generated catch block
			   e1.printStackTrace();
		   } catch (InvocationTargetException e1) {
			   // TODO Auto-generated catch blockPotentialJList
			   e1.printStackTrace();
		   } catch (NoSuchMethodException e1) {
			   // TODO Auto-generated catch block
			   e1.printStackTrace();
		   }
        }
	}

	

	
	
	
	
}
