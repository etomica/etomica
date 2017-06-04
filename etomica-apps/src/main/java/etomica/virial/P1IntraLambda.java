/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.IPotentialAtomic;
import etomica.potential.P1IntraMolecular;
import etomica.space.Space;

public class P1IntraLambda implements IPotentialAtomic, P1IntraMolecular {
	protected double lambda = -1, u0 = 0;
	protected P1IntraMolecular p1;
	public P1IntraLambda(Space space, double lambda, P1IntraMolecular pot1, double u00) {
		p1 = pot1;
		this.lambda = lambda;		        
		u0 = u00;
    }

	public double getRange() {
		return 0;
	}

	public void setBox(Box box) {
		
	}

	public int nBody() {
		return 1;
	}

	public double energy(IAtomList atoms) {
		if (lambda == -1) throw new RuntimeException("lambda needs to be set first");
		return lambda*(((IPotentialAtomic)p1).energy(atoms) - u0);
	}
	
	public double du(double r) {	
		return lambda*p1.du(r);
	}

	public double d2u(double r) {		
		return lambda*p1.d2u(r);
	}
	
	public double u(double r) {
		return lambda*(p1.u(r) - u0);
	}
	
}
