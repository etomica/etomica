/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.atom.IMoleculeList;
import etomica.api.IPotential;

public class MayerFunctionSum implements MayerFunction {

	public MayerFunctionSum(MayerFunction[] functions, double []coefficients) {
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

	public double f(IMoleculeList pair, double r2, double beta) {
		double sum = 0;
		for (int i=0;i<functions.length;i++){
			sum += coefficients[i]*functions[i].f(pair, r2, beta);
		}
		return sum;
	}

	protected final MayerFunction[] functions;
	protected final double []coefficients;

}
