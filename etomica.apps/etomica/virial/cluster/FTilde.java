/*
 * Created on Mar 12, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial.cluster;

import etomica.space.CoordinatePair;
import etomica.virial.MayerFunction;


public class FTilde implements MayerFunction {
	private final MayerFunction fWrapped;
	public FTilde(MayerFunction f) {
		fWrapped = f;
	}
	public double f(CoordinatePair cPair, double beta) {
		return fWrapped.f(cPair,beta) + 1.0;
	}
	public String toString() {return "f~  ";}
}