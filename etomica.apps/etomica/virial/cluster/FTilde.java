package etomica.virial.cluster;

import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IPotential;

import etomica.virial.MayerFunction;


public class FTilde implements MayerFunction, java.io.Serializable {
	private final MayerFunction fWrapped;
	public FTilde(MayerFunction f) {
		fWrapped = f;
	}
	public double f(IAtomSet aPair, double beta) {
		return fWrapped.f(aPair,beta) + 1.0;
	}
	public IPotential getPotential() {
	    return fWrapped.getPotential();
	}

    public void setBox(IBox newBox) {
        fWrapped.setBox(newBox);
	}

	public String toString() {return "f~  ";}
}
