package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotential;
import etomica.potential.P2HardAssociationCone;
import etomica.space.ISpace;

public class MayerFunctionProductGeneral implements MayerFunction {

	public MayerFunctionProductGeneral(ISpace space, MayerFunction[]functions, double []coefficients) {
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
		double prod = 1;
		for (int i=0;i<functions.length;i++){
			
			prod *= coefficients[i]*functions[i].f(pair, r2, beta);
			//if (P2HardAssociationCone.FLAG)System.out.println (i+" " + functions[i].f(pair, r2, beta));
			}
		//if (P2HardAssociationCone.FLAG)System.exit(0);
		return prod;
	}

}
