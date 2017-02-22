/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;

import etomica.chem.elements.ElementSimple;

public class SimpleElementForSimilarSpecies {
	
	

	public ElementSimple CH3element;
    public ElementSimple CH2element;
    public ElementSimple TSelement;
    public ElementSimple Aelement;
    public ElementSimple Oelement;
    
    
    public ElementSimple getOelement() {
		return Oelement;
	}

	public void setOelement(ElementSimple oelement) {
		Oelement = oelement;
	}

	private SimpleElementForSimilarSpecies(){
    	this.CH3element = new ElementSimple("CH3",15);
    	this.CH2element = new ElementSimple("CH2",14);
    	this.TSelement = new ElementSimple("TS");
    	this.Aelement = new ElementSimple("A");
    	this.Oelement = new ElementSimple("O");
    }
    
    /**
     * SingletonHolder is loaded on the first execution of Singleton.getInstance() 
     * or the first access to SingletonHolder.INSTANCE, not before.
     */
     private static class SimpleElementForAlkaneSpeciesHolder { 
             public static final SimpleElementForSimilarSpecies instance = new SimpleElementForSimilarSpecies();
     }
     
     public ElementSimple getAelement() {
		return Aelement;
	}

	public void setAelement(ElementSimple aelement) {
		Aelement = aelement;
	}

	public static SimpleElementForSimilarSpecies getInstance() {
         return SimpleElementForAlkaneSpeciesHolder.instance;
     }
     
     public ElementSimple getCH3element() {
 		return CH3element;
 	}

 	public void setCH3element(ElementSimple cH3element) {
 		CH3element = cH3element;
 	}

 	public ElementSimple getCH2element() {
 		return CH2element;
 	}

 	public void setCH2element(ElementSimple cH2element) {
 		CH2element = cH2element;
 	}

	public ElementSimple getTSelement() {
		return TSelement;
	}

	public void setTSelement(ElementSimple tSelement) {
		TSelement = tSelement;
	}
 	
 	
     
}
