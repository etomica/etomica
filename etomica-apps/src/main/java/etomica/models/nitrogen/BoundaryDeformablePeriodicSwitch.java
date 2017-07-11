/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.space.Vector;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;

public class BoundaryDeformablePeriodicSwitch extends BoundaryDeformablePeriodic{
	
	public BoundaryDeformablePeriodicSwitch(Space _space, Vector[] vex) {
		super(_space, vex);
		temp = space.makeVector();
		
	}
	public void nearestImage(Vector dr) {
		if(doPBC){
			super.nearestImage(dr);
		} 
	}

	public Vector centralImage(Vector r) {
		
		if(doPBC){
			return super.centralImage(r);
		}
		temp.E(0.0);
		return temp;
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
	protected Vector temp;
	protected boolean doPBC = true;
}
