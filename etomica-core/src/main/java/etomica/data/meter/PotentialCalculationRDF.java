package etomica.data.meter;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.api.IVectorMutable;
import etomica.data.DataSourceUniform;
import etomica.data.DataSourceUniform.LimitType;
import etomica.potential.PotentialCalculation;
import etomica.space.ISpace;
import etomica.units.Length;

/**
 * PotentialCalculation that simply does the work of collecting an RDF.
 *
 * @author Andrew Schultz
 */

public class PotentialCalculationRDF implements PotentialCalculation {

    protected final IVectorMutable dr;
    protected long[] gSum;
    protected IBoundary boundary;
    protected final DataSourceUniform xDataSource;
    protected double xMax;

    
    public PotentialCalculationRDF(ISpace space, IBox box) {
        dr = space.makeVector();
        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
        
        gSum = new long[xDataSource.getData().getLength()];
        
        this.boundary = box.getBoundary();
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }
    
    public void setBox(IBox box) {
        boundary = box.getBoundary();
    }

    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        xMax = xDataSource.getXMax();
        gSum = new long[xDataSource.getData().getLength()];
    }

    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        IAtom atom0 = atoms.getAtom(0);
        IAtom atom1 = atoms.getAtom(1);
        dr.Ev1Mv2(atom0.getPosition(), atom1.getPosition());

        double xMaxSquared = xMax*xMax;
        boundary.nearestImage(dr);
        double r2 = dr.squared();       //compute pair separation
        if(r2 < xMaxSquared) {
            int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
            gSum[index]++;                        //add once for each atom
        }
    }

    public long[] getGSum() {
        return gSum;
    }

}
