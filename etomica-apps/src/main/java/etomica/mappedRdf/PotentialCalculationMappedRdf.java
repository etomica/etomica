package etomica.mappedRdf;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceUniform;
import etomica.potential.IPotentialAtomic;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;

public class PotentialCalculationMappedRdf implements PotentialCalculation {
    protected final Vector dr;
    protected double[] gSum;
    protected double[] gSum2;
    protected double[] newgSum;
    protected double[] newestgSum;
    protected double[] thirdterm;
    protected double rcforHandfinmap;
    protected Boundary boundary;
    protected final DataSourceUniform xDataSource;
    protected double xMax;
    protected double vol;
    protected double q;
    protected double qp;
    protected double sum, R;
    protected double x0, vCut;
    protected double vShift;
    protected double beta;
    protected double c1;
    protected final int nbins;
    protected final Space space;
    protected final double[] cumint;
    protected final AtomLeafAgentManager<Vector> forceManager;
    protected Potential2SoftSpherical p2;


    public PotentialCalculationMappedRdf(double rcforHandfinmap, Space space, Box box, int nbins, AtomLeafAgentManager<Vector> forceManager) {
        dr = space.makeVector();
        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(DataSourceUniform.LimitType.HALF_STEP);
        xDataSource.setTypeMin(DataSourceUniform.LimitType.HALF_STEP);
        newgSum=new double[xDataSource.getData().getLength()];
        newestgSum=new double[xDataSource.getData().getLength()];

        gSum = new double[xDataSource.getData().getLength()];
        gSum2 = new double[xDataSource.getData().getLength()];

        thirdterm = new double[xDataSource.getData().getLength()];
        this.rcforHandfinmap = rcforHandfinmap;

        this.boundary = box.getBoundary();
        this.nbins = nbins;
        this.space = space;

        cumint = new double[nbins + 1];
        vol = box.getBoundary().volume();

        this.forceManager = forceManager;
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double getX0() {
        return x0;
    }

    public double getq() {
        return q;
    }

    public void setVCut(double newVCut) {
        vCut = newVCut;
    }

    /**
     * Sets volume to an arbitrary value (instead of the box volume).  This is
     * used for 1/(1+q/V) terms.
     */
    public void setVolume(double newVol) {
        vol = newVol;
    }

    public void setPotential(Potential2SoftSpherical p) {
        p2 = p;
    }

    public void setTemperature(double T) {
        beta = 1 / T;
        double rc = p2.getRange();
        x0 = rc;
        System.out.println("x0 " + x0);
        q = 0;

        if (vCut == 0) vCut = x0;
        vShift = -p2.u(vCut * vCut);
        c1 = Math.log(rc + 1) / nbins;
        int D = space.D();
        for (int i = 1; i <= nbins; i++) {
            double r = Math.exp(c1 * i) - 1;
            double r2 = r * r;
            if (r >= p2.getRange()) {
                if (i < nbins) throw new RuntimeException("oops " + i + " " + nbins + " " + r + " " + p2.getRange());
                r = Math.exp(c1 * i * 0.9999) - 1;
                r2 = r * r;
            }
            double u = p2.u(r2);
            double v = calcV(r, u);
            double evm1 = Math.exp(-beta * v);
            q += (D == 2 ? r : r2) * (evm1 - 1) * c1 * (r + 1);

        }

        q *= (D == 2 ? 2 : 4) * Math.PI;
        q += vol;

        for (int i = 1; i <= nbins; i++) {
            double r = Math.exp(c1 * i) - 1;
            double r2 = r * r;
            if (r >= p2.getRange()) {
                if (i < nbins) throw new RuntimeException("oops " + i + " " + nbins + " " + r + " " + p2.getRange());
                r = Math.exp(c1 * i * 0.9999) - 1;
                r2 = r * r;
            }
            double u = p2.u(r2);
            double v = calcV(r, u);
            double evm1 = Math.exp(-beta * v);
            cumint[i] = cumint[i - 1] + (D == 2 ? r : r2) * (evm1) * c1 * (r + 1);

        }

    }

    protected double calcV(double r, double u) {

        if (r > vCut)
            return 0;
        return 0; //u+vShift;
    }

    public double[] gR() {
        double[] gR = new double[xDataSource.getData().getLength()];
        for (int k = 0; k < xDataSource.getData().getLength(); k++) {
            double R = xDataSource.getData().getValue(k);
            double uR = p2.u(R);
            double vR = calcV(R, uR);
            //  System.out.println(vR*beta);
            double evmR = Math.exp(-beta * vR);
            gR[k] = evmR * vol / q;

        }
        return gR;
    }

    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        xMax = xDataSource.getXMax();
        gSum = new double[xDataSource.getData().getLength()];
        newgSum = new double[xDataSource.getData().getLength()];
        newestgSum = new double[xDataSource.getData().getLength()];

        gSum2 = new double[xDataSource.getData().getLength()];
        thirdterm= new double[xDataSource.getData().getLength()];
    }

    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof Potential2SoftSpherical)) return;
        IAtom atom0 = atoms.get(0);
        IAtom atom1 = atoms.get(1);
        dr.Ev1Mv2(atom0.getPosition(), atom1.getPosition());

        double xMaxSquared = xMax * xMax;
        boundary.nearestImage(dr);
        double r2 = dr.squared();       //compute pair separation
        double r = Math.sqrt(r2);
        if (r > rcforHandfinmap) return;//<l/2

        if (r2 < xMaxSquared) {
            for (int k = 0; k < xDataSource.getData().getLength(); k++) {
                double R = xDataSource.getData().getValue(k);
                Vector fi = forceManager.getAgent(atom0);
                Vector fj = forceManager.getAgent(atom1);
                     double xu = R < r ? 1 : 0; //calcXu(r, u, R);
               //      gSum[k] -= ((xu/(4 * Math.PI)))* wp*beta;               //add once for each atom
                     //replaced fi.dot(dr)-fj.dot(dr) terms with fi.dot(dr)
                    gSum[k] =gSum[k] - ((xu/(4 * Math.PI*r * r * r))-(1/(3*vol)))* (fi.dot(dr)-fj.dot(dr)  )*beta/vol;               //add once for each atom
                    //   gSum[k] -= ((-r*r*r/(3*vol)))* wp*beta;               //add once for each atom
                     gSum2[k] = gSum2[k] - ((xu/(4 * Math.PI*r * r * r)) )* (fi.dot(dr)-fj.dot(dr)  )*beta/vol;                //add once for each atom
                     thirdterm[k] +=  (beta*(fi.dot(dr)-fj.dot(dr)  )/(3*vol*(4*Math.PI*rcforHandfinmap*rcforHandfinmap*rcforHandfinmap/3)));               //add once for each atom

                newgSum[k]=newgSum[k]-((beta/vol)*-1*fi.dot(dr)*((xu/(4*Math.PI*r * r * r))-(1/(3*vol))));

                     newestgSum[k]=newestgSum[k]+(beta*-1*fi.dot(dr)*(1-xu)/(vol*4*Math.PI*r * r * r));
//
            }
        }
     }

    public double[] getGSum() {
        return gSum;
    }

}
