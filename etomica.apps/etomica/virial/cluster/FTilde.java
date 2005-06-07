package etomica.virial.cluster;

import etomica.AtomSet;
import etomica.virial.MayerFunction;


public class FTilde implements MayerFunction {
	private final MayerFunction fWrapped;
	public FTilde(MayerFunction f) {
		fWrapped = f;
	}
	public double f(AtomSet aPair, double beta) {
		return fWrapped.f(aPair,beta) + 1.0;
	}
	public String toString() {return "f~  ";}
}