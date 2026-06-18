/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;

public class ConformationBenzene implements IConformation {

	public ConformationBenzene(double bondL){
		this.bondL = bondL;
	}

	@Override
	public void initializePositions(IAtomList atomList) {
		double radius = bondL;
		for (int i=0; i<6; i++) {
			IAtom a = atomList.get(i);
			Vector p = a.getPosition();
			p.setX(0, radius*Math.cos(Math.PI*i/3));
			p.setX(1, radius*Math.sin(Math.PI*i/3));
			p.setX(2, 0);
		}
	}

	protected double bondL;
}


