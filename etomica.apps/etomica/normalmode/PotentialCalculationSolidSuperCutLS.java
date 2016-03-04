package etomica.normalmode;

import etomica.api.IAtomList;
import etomica.api.IPotentialAtomic;
import etomica.liquidLJ.Potential2SoftSphericalLSMultiLat;
import etomica.space.ISpace;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationSolidSuperCutLS extends PotentialCalculationSolidSuperCut {
        
    public PotentialCalculationSolidSuperCutLS(ISpace space, CoordinateDefinition coordinateDefinition, double[] cutoffs) {
        super(space, coordinateDefinition, cutoffs);
    }
    
    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof Potential2SoftSphericalLSMultiLat)) return;
        Potential2SoftSphericalLSMultiLat potentialSoft = (Potential2SoftSphericalLSMultiLat)potential;
        double[][] ijSums = potentialSoft.energyVirialCut(atoms);

        int n = ijSums[0].length;
        if (n != energySum.length) {
            energySum = new double[n];
            virialSum = new double[n];
            sum1 = new double[n];
            dadbSum = new double[n];
        }
        for (int i=0; i<n; i++) {
            energySum[i] += ijSums[0][i];
            virialSum[i] += ijSums[1][i];
            sum1[i] += ijSums[2][i]*fac1;
            dadbSum[i] += ijSums[3][i];
        }
        
    }
}
