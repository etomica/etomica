/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.components;

public class ExceptionDoNotRemoveSpecies extends Exception{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public int indexSelectedSpecies;
	
	public ExceptionDoNotRemoveSpecies(int index){
		indexSelectedSpecies = index;
	}

	public int getIndexSelectedSpecies() {
		return indexSelectedSpecies;
	}

	public void setIndexSelectedSpecies(int indexSelectedSpecies) {
		this.indexSelectedSpecies = indexSelectedSpecies;
	}
	
	
}
