package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotential;
import etomica.space.ISpace;

public class MayerFunctionSumGeneral implements MayerFunction {

	public MayerFunctionSumGeneral(ISpace space, MayerFunction[]functions, double []coefficients) {
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
