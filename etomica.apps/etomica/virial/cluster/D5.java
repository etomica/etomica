package etomica.virial.cluster;

import etomica.virial.Cluster;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * The virial cluster of 4 points joined as a simple ring, with one bond
 * joining an opposite pair.
 */
public final class D5 extends Cluster {
	public D5(MayerFunction f) {
		super(4, -3./4., new BondGroup(f, Standard.D5));
	}
}
