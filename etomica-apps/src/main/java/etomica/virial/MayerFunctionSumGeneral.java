/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IPotential;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;

public class MayerFunctionSumGeneral implements MayerFunction {

	public MayerFunctionSumGeneral(Space space, MayerFunction[]functions, double []coefficients) {
		this.functions = functions;
		this.coefficients = coefficients;
	}

	public IPotential getPotential() {
		return null;
	}

	public void setBox(Box box) {
		for (int i=0;i<functions.length;i++){
			functions[i].setBox(box);
			}

	}

	
	protected final MayerFunction[]functions;
	protected final double []coefficients;
	public double f(IMoleculeList pair, double r2, double beta) {
		double sum = 0;
		for (int i=0;i<functions.length;i++){
			
			sum += coefficients[i]*functions[i].f(pair, r2, beta);
			}
		return sum;
	}

}
