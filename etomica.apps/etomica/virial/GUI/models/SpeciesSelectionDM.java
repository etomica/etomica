package etomica.virial.GUI.models;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;

public class SpeciesSelectionDM {
	
	private static final String[] intialDisplayList = {"---No species selected---","Press \'Add\' to select a species"};
	private static final String[] intialPotentialDisplayList = {"---No potentials selected---","Press \'Add\' to select a potential"};
	
	private static final String[] species_List = {"LJ","CO2","Methanol","Ethanol","Methane","Ethane","Propane","Higher n-Alkanes","Naphthalene","Water"};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesLJ_PotentialsList = {CreateSpeciesDM_LJ_LJ.class,CreateSpeciesDM_LJ_LJQ.class,CreateSpeciesDM_LJ_2CLJQ.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesCO2_PotentialsList = {CreateSpeciesDM_CO2_2CLJQ.class,CreateSpeciesDM_CO2_EMP2.class,CreateSpeciesDM_CO2_Trappe.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesMethane_PotentialList = {CreateSpeciesDM_Alkane_TRAPPE.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesEthane_PotentialList = {CreateSpeciesDM_Alkane_TRAPPE.class, CreateSpeciesDM_Alkane_SKS.class, CreateSpeciesDM_Ethane_2CLJQ.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesPropane_PotentialList = {CreateSpeciesDM_Alkane_TRAPPE.class,CreateSpeciesDM_Alkane_SKS.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesNAlkane_PotentialList = {CreateSpeciesDM_Alkane_TRAPPE.class, CreateSpeciesDM_Alkane_SKS.class};
	
	@SuppressWarnings("rawtypes")
	private static final Class[] speciesH2O_PotentialList = {CreateSpeciesDM_H2O_SPCE.class};
	
	private String speciesIdentifier;
	private int[] speciesSelectedIndex;
	private int[] potentialsSelectedIndex;
	private String[] potentialsList;
	
	public SpeciesSelectionDM(){
		reset();
	}
	
	public void reset() {
		speciesIdentifier = null;
		speciesSelectedIndex = new int[2];
		potentialsSelectedIndex = new int[2];
	}
	
	public String getSpeciesIdentifier() {
		return speciesIdentifier;
	}

	public void setSpeciesIdentifier(String speciesIdentifier) {
		this.speciesIdentifier = speciesIdentifier;
	}

	public String[] getPotentialsList() {
		return potentialsList;
	}

	public int[] getSpeciesSelectedIndex() {
		return speciesSelectedIndex;
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
	
	public void setPotentialList(Class[] potentialList){
		this.potentialsList = new String[potentialList.length];
        for (int i=0; i<potentialsList.length; i++) {
		   try {
			   try {
				   potentialsList[i] = (String)potentialList[i].getMethod("getCustomName", new Class[0]).invoke(potentialList[i].getConstructor().newInstance(new Object[0]),new Object[0]);
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

	public static Class[] getSpeciesh2oPotentiallist() {
		return speciesH2O_PotentialList;
	}

	
	
	
	
}
