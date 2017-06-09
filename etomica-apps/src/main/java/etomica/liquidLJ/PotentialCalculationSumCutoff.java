package etomica.liquidLJ;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Evaluates the energy summed over all iterated atoms. Each call to doCalculate
 * accumulates the energy of the given potential applied to the atoms produced
 * by the given iterator. Each such call accumulates a single sum, which can
 * be accessed via the getSum() method.  Sum is re-zeroed only upon a call to the
 * zeroSum() method. 
 *
 * @author David Kofke
 */
public class PotentialCalculationSumCutoff implements PotentialCalculation {

    public PotentialCalculationSumCutoff(Space space, double[] cutoffs) {
        dr = space.makeVector();
        r2Cuts = new double[cutoffs.length];
        for (int i=0; i<cutoffs.length; i++) {
            r2Cuts[i] = cutoffs[i]*cutoffs[i];
        }
        uSums = new double[cutoffs.length];
        vSums = new double[cutoffs.length];
    }
    
    public void setBox(Box box) {
        this.box = box;
        boundary = box.getBoundary();
    }

    /**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 > r2Cuts[r2Cuts.length-1]) return;
        double u = ((Potential2SoftSpherical)potential).u(r2);
        double v = ((Potential2SoftSpherical)potential).du(r2);
        for (int i=uSums.length-1; i>=0; i--) {
            if (r2 > r2Cuts[i]) break;
            uSums[i] += u;
            vSums[i] += v;
        }
	}
	
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 */
	public void zeroSums() {
	    for (int i=0; i<vSums.length; i++) {
	        uSums[i] = vSums[i] = 0.0;
	    }
		boundary = box.getBoundary();
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double[] getUSums() {
        return uSums;
    }
	
	/**
     * Returns the current value of the energy sum.
     */
    public double[] getVSums() {
        return vSums;
    }

	
	protected double[] uSums, vSums, r2Cuts;
	protected final Vector dr;
	protected Box box;
	protected Boundary boundary;
}
