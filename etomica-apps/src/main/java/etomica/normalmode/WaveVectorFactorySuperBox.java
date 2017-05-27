/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.space.Vector;
import etomica.box.Box;
import etomica.lattice.crystal.Primitive;
import etomica.space.Space;

public class WaveVectorFactorySuperBox extends WaveVectorFactorySimple {

	/**
	 * Specifically for FCC 32 system size
	 * 
	 * Instead of making wavevectors for 216 cells; this class will
	 * 	reduce the size of box's boundary by one-third in order 
	 * 	to make wave vectors of 8 cells
	 * 
	 * and then the size of box's boundary is restored to the initial
	 * 	condition. 
	 * 
	 * @author Tai Boon Tan
	 */
	private static final long serialVersionUID = 1L;

	public WaveVectorFactorySuperBox(Primitive primitive, Space _space) {
		super(primitive, _space);
		// TODO Auto-generated constructor stub
	}
	
	public void makeWaveVectors(Box box){
		
		Vector boxDimension = space.makeVector();
		boxDimension.E(box.getBoundary().getBoxSize());
		double fraction = 1.0/2;
		boxDimension.TE(fraction);
		box.getBoundary().setBoxSize(boxDimension);
		
		super.makeWaveVectors(box);
		
		boxDimension.TE(2.0);
		box.getBoundary().setBoxSize(boxDimension);
		
	}
	
}
