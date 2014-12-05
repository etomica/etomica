/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;

import java.util.ArrayList;

import etomica.virial.GUI.models.MolecularModelTRAPPE_SpeciesAlkane;
import etomica.virial.GUI.models.MolecularModel2CLJQ_SpeciesCO2;
import etomica.virial.GUI.models.MolecularModelEMP2_SpeciesCO2;
import etomica.virial.GUI.models.MolecularModelTrappe_SpeciesCO2;
import etomica.virial.GUI.models.MolecularModel2CLJQ_SpeciesLJ;
import etomica.virial.GUI.models.MolecularModelLJ_SpeciesLJ;
import etomica.virial.GUI.models.MolecularModelLJQ_SpeciesLJ;




public class ArrayMoleculeDisplayList {
	
	


	private static int Index=0;
	private int id;
	private Object[] SpeciesL;
	
	public ArrayMoleculeDisplayList(){
		id = Index;
		SpeciesL = new Object[8];
	}
	
	public void addSpecies(Object object){
		SpeciesL[id] = object;
		id++;
	}
	
	
	public void removeSpecies(){
		SpeciesL[id - 1] = null;
		if(id > 0){
			id--;
		}
		else{
			id = 0;
		}
	}
	public void removeSpeciesAtIndex(int index){
		SpeciesL[index] = null;
		if(id > 0){
			id--;
		}
		else{
			id = 0;
		}
	}
	
	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}
	
	public Object getObject(int Index){
		return SpeciesL[Index];
	}
	
	public ArrayList<String> displayList(){
		ArrayList<String> DisplayArray = new ArrayList<String>();
		
		for (int i=0;i<id;i++){
			if(SpeciesL[i] != null){
				if(SpeciesL[i] instanceof MolecularModelLJ_SpeciesLJ){
					DisplayArray.add("LJ-Spherical2Body");
				}
				else if(SpeciesL[i] instanceof MolecularModelLJQ_SpeciesLJ){
					DisplayArray.add("LJ-Spherical2BodyWithQ");
				}
				else if(SpeciesL[i] instanceof MolecularModel2CLJQ_SpeciesLJ){
					DisplayArray.add("2CenteredLJWithQ");
				}
				
				else if(SpeciesL[i] instanceof MolecularModel2CLJQ_SpeciesCO2){
					DisplayArray.add("CO2-2CenteredLJWithQ");
				}
				else if(SpeciesL[i] instanceof MolecularModelEMP2_SpeciesCO2){
					DisplayArray.add("CO2-EPM2");
				}
				else if(SpeciesL[i] instanceof MolecularModelTrappe_SpeciesCO2){
					DisplayArray.add("CO2-TRAPPE");
				}
				else if(SpeciesL[i] instanceof MolecularModelTRAPPE_SpeciesAlkane){
					DisplayArray.add("n-Alkane-TRAPPE");
				}
				
			}
		}
		
		return DisplayArray;
	}
}
