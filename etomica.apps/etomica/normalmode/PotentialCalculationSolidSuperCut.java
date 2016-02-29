package etomica.normalmode;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.space.ISpace;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationSolidSuperCut implements PotentialCalculation {
        
    protected final CoordinateDefinition coordinateDefinition;
    protected final IVectorMutable drSite0, drSite1, drA, dr, drB;
    protected final ISpace space;
    protected double[] sum1, virialSum;
    protected double[] energySum, dadbSum;
    protected double fac1;
    protected IBox box;
    protected IBoundary boundary;
    protected final double[] r2Cut;
    
    public PotentialCalculationSolidSuperCut(ISpace space, CoordinateDefinition coordinateDefinition, double[] cutoffs) {
        this.coordinateDefinition = coordinateDefinition;
        this.space = space;
        this.r2Cut = new double[cutoffs.length];
        for (int i=0; i<r2Cut.length; i++) {
            r2Cut[i] = cutoffs[i]*cutoffs[i];
        }
        drSite0 = space.makeVector();
        drSite1 = space.makeVector();
        drA = space.makeVector();
        dr = space.makeVector();
        drB = space.makeVector();
        setBox(coordinateDefinition.getBox());
        sum1 = new double[r2Cut.length];
        virialSum = new double[r2Cut.length];
        energySum = new double[r2Cut.length];
        dadbSum = new double[r2Cut.length];
    }
    
    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        IAtom atom0 = atoms.getAtom(0);
        IAtom atom1 = atoms.getAtom(1);
        IVector site0 = coordinateDefinition.getLatticePosition(atom0);
        IVector site1 = coordinateDefinition.getLatticePosition(atom1);

        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        Potential2SoftSpherical potentialSoft = (Potential2SoftSpherical)potential;

        drSite0.Ev1Mv2(atom0.getPosition(), site0);
        drSite1.Ev1Mv2(atom1.getPosition(), site1);
        drA.Ev1Mv2(drSite1, drSite0);

        drB.Ev1Mv2(site1, site0);
        boundary.nearestImage(drB);
        double r2site = drB.squared();
        
        double u = potentialSoft.u(r2);
        double du = potentialSoft.du(r2);
        double p1 = -du/r2*dr.dot(drB)*fac1;
        double dadb = -du/r2*dr.dot(drA);

        int n=r2Cut.length;
        for (int i=n-1; i>=0; i--) {
            if (r2Cut[i] < r2site) break;
            energySum[i] += u;
            sum1[i] += p1;
            virialSum[i] += du;
            dadbSum[i] += dadb;
        }
        
    }
    
    public double[] getPressure1() {
        return sum1;
    }
    
    public double[] getVirialSum() {
        return virialSum;
    }

    public double[] getEnergySum() {
      return energySum;
    }
    
    /**
     * Returns sum of Fi dot ri
     */
    public double[] getDADBSum() {
      return dadbSum;
    }
    
    public void zeroSum() {
        int n = r2Cut.length;
        for (int i=0; i<n; i++) {
            sum1[i] = virialSum[i] = energySum[i] = dadbSum[i] = 0;
        }
        boundary = box.getBoundary();
        int D = boundary.getBoxSize().getD();
        fac1 = 1.0/(D*boundary.volume());
    }
    
    public void setBox(IBox box) {
        this.box = box;
        boundary = box.getBoundary();
    }
}
