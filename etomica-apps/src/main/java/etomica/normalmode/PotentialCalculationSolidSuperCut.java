package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.api.IBoundary;
import etomica.box.Box;
import etomica.api.IPotentialAtomic;
import etomica.space.Vector;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationSolidSuperCut implements PotentialCalculation {
        
    protected final CoordinateDefinition coordinateDefinition;
    protected final Vector drSite0, drSite1, drA, dr, drB;
    protected final Space space;
    protected Vector[] pSumXYZ1, pSumXYZ2;
    protected final Vector pTmp1, pTmp2;
    protected double[] sum1, virialSum;
    protected double[] energySum, dadbSum;
    protected double fac1;
    protected double[] fac2;
    protected Box box;
    protected IBoundary boundary;
    protected final double[] r2Cut;
    protected double temperature;
    protected double[] pHarmonic;
    
    public PotentialCalculationSolidSuperCut(Space space, CoordinateDefinition coordinateDefinition, double[] cutoffs) {
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
        pSumXYZ1 = new Vector[r2Cut.length];
        pSumXYZ2 = new Vector[r2Cut.length];
        for (int i=0; i<pSumXYZ1.length; i++) {
            pSumXYZ1[i] = space.makeVector();
            pSumXYZ2[i] = space.makeVector();
        }
        pTmp1 = space.makeVector();
        pTmp2 = space.makeVector();
    }
    
    public void setHarmonicPressure(double[] pHarmonic, double temperature) {
        this.pHarmonic = pHarmonic;
        this.temperature = temperature;
    }
    
    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        IAtom atom0 = atoms.getAtom(0);
        IAtom atom1 = atoms.getAtom(1);
        Vector site0 = coordinateDefinition.getLatticePosition(atom0);
        Vector site1 = coordinateDefinition.getLatticePosition(atom1);

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
        pTmp1.Ea1Tv1(-du/r2*fac1, dr);
        pTmp1.TE(drB);
        double dadb = -du/r2*dr.dot(drA);
        pTmp2.Ea1Tv1(-du/r2, dr);
        pTmp2.TE(drA);

        int n=r2Cut.length;
        for (int i=n-1; i>=0; i--) {
            if (r2Cut[i] < r2site) break;
            energySum[i] += u;
            sum1[i] += p1;
            virialSum[i] += du;
            dadbSum[i] += dadb;
            pSumXYZ1[i].PE(pTmp1);
            pSumXYZ2[i].PE(pTmp2);
        }
        
    }
    
    public double[] getPressure1() {
        return sum1;
    }
    
    public Vector[] getPressure1XYZ() {
        return pSumXYZ1;
    }
    
    public Vector[] getDADBXYZ() {
        return pSumXYZ2;
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
            pSumXYZ1[i].E(0);
            pSumXYZ2[i].E(0);
        }
        boundary = box.getBoundary();
        int D = boundary.getBoxSize().getD();
        double vol = boundary.volume();
        fac1 = 1.0/(D*vol);
    }
    
    public void setBox(Box box) {
        this.box = box;
        boundary = box.getBoundary();
    }
}
