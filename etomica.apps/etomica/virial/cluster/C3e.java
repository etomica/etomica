package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.MayerFunction;
import etomica.virial.PairSet;

/**
 * @author kofke
 *
 * e-bond cluster: e12 e13 e23 - (e12+e13+e23) + 2
 */
/* History
 * 08/21/03 (DAK) removed reset() calls for pairs in value method.  Reset should
 * be done by class invoking this method.
 */
public class C3e extends Cluster {

	public C3e(MayerFunction f) {
		super(3, -1./3.,new BondGroup(f, Standard.C3));
	}


	
	/**
	 * @see etomica.virial.Cluster#value(etomica.virial.PairSet, double)
	 */
	public double value(PairSet pairs, double beta) {
		double e12 = bondArray[0][1].f(pairs.getPair(0,1),beta) + 1.0;
		double e13 = bondArray[0][2].f(pairs.getPair(0,2),beta) + 1.0;
		double e23 = bondArray[1][2].f(pairs.getPair(1,2),beta) + 1.0;
		return e12*e13*e23 - 1;// - (e12 + e13 + e23) + 2;
	}

}
