package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.space.ISpace;

public class MayerFunctionSum extends MayerFunctionSpherical {

	public MayerFunctionSum(ISpace space, MayerFunctionSpherical[]functions, double []coefficients) {
		super(space);
		this.functions = functions;
		this.coefficients = coefficients;
	}

	public IPotential getPotential() {
		return null;
	}

	public void setBox(IBox box) {
		for (int i=0;i<functions.length;i++){
			functions[i].setBox(box);
			}

	}

	public double f(double r2, double beta) {
		double sum = 0;
		for (int i=0;i<functions.length;i++){
			
			sum += coefficients[i]*functions[i].f(r2, beta);
			}
		return sum;
	}
	protected final MayerFunctionSpherical[]functions;
	protected final double []coefficients;

}
