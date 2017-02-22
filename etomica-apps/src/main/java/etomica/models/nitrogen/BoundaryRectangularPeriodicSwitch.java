/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;

public class BoundaryRectangularPeriodicSwitch extends BoundaryRectangularPeriodic{
	
	public BoundaryRectangularPeriodicSwitch(ISpace _space) {
		super(_space);
		
	}
	public void nearestImage(IVectorMutable dr) {
		if(doPBC){
			super.nearestImage(dr);
		} 
	}

	public IVector centralImage(IVector r) {
		
		if(doPBC){
			return super.centralImage(r);
		}
		tempImage.E(0.0);
		return tempImage;
	}

	public boolean getPeriodicity(int d) {
		
		return doPBC;
	}

	public boolean isDoPBC() {
		return doPBC;
	}

	public void setDoPBC(boolean doPBC) {
		this.doPBC = doPBC;
	}

	private static final long serialVersionUID = 1L;
	protected boolean doPBC = true;
}
