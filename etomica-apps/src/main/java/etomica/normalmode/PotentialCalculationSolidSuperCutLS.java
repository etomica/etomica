package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.api.IPotentialAtomic;
import etomica.space.Vector;
import etomica.liquidLJ.Potential2SoftSphericalLSMultiLat;
import etomica.liquidLJ.Potential2SoftSphericalLSMultiLat.ReturnValue;
import etomica.space.Space;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationSolidSuperCutLS extends PotentialCalculationSolidSuperCut {
        
    public PotentialCalculationSolidSuperCutLS(Space space, CoordinateDefinition coordinateDefinition, double[] cutoffs) {
        super(space, coordinateDefinition, cutoffs);
    }
    
    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof Potential2SoftSphericalLSMultiLat)) return;
        Potential2SoftSphericalLSMultiLat potentialSoft = (Potential2SoftSphericalLSMultiLat)potential;
        ReturnValue rv = potentialSoft.energyVirialCut(atoms);
        int n = rv.energySum.length;
        if (n != energySum.length) {
            energySum = new double[n];
            virialSum = new double[n];
            sum1 = new double[n];
            dadbSum = new double[n];
            pSumXYZ1 = new Vector[n];
            pSumXYZ2 = new Vector[n];
            for (int i=0; i<n; i++) {
                pSumXYZ1[i] = space.makeVector();
                pSumXYZ2[i] = space.makeVector();
            }
        }
        for (int i=0; i<n; i++) {
            energySum[i] += rv.energySum[i];
            virialSum[i] += rv.virialSum[i];
            sum1[i] += rv.sum1[i]*fac1;
            dadbSum[i] += rv.dadbSum[i];
            pSumXYZ1[i].PEa1Tv1(fac1, rv.pSumXYZ1[i]);
            pSumXYZ2[i].PE(rv.pSumXYZ2[i]);
        }
    }
}
