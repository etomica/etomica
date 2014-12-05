/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

public class ModelAlkaneLengthSelection {
	
	private int[] alkaneLength;
	
	
	public ModelAlkaneLengthSelection(){
		reset();
	}

	public void reset() {
		// TODO Auto-generated method stub
		alkaneLength = new int[10];
		for(int i=0;i<10;i++){
			alkaneLength[i]=4;
		}
		
		
	}

	public int getAlkaneLength(int index) {
		return alkaneLength[index];
	}

	public void setAlkaneLength(int alkane1Length,int index) {
		this.alkaneLength[index] = alkane1Length;
	}

	

}
