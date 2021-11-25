package etomica.mappedRdf;

import etomica.box.Box;
import etomica.data.DataSourceUniform;
import etomica.potential.Potential2Soft;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.units.dimensions.Length;

public class PotentialCallbackMappedRdf implements PotentialCallback {

    protected final PotentialCompute potentialCompute;
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
    protected double sum;
    protected double x0, vCut;
    protected double vShift;
    protected double beta;
    protected double c1;
    protected final int nbins;
    protected final double[] cumint;
    protected Potential2Soft p2;


    public PotentialCallbackMappedRdf(double rcforHandfinmap, Box box, int nbins, PotentialCompute potentialCompute) {
        this.potentialCompute = potentialCompute;
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

        cumint = new double[nbins + 1];
        vol = box.getBoundary().volume();
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

    public void setPotential(Potential2Soft p) {
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
        int D = boundary.getBoxSize().getD();
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

    @Override
    public void pairCompute(int i, int j, Vector dr, double[] u012) {

        double xMaxSquared = xMax * xMax;
        double r2 = dr.squared();       //compute pair separation
        double r = Math.sqrt(r2);
        if (r > rcforHandfinmap) return;//<l/2

        if (r2 < xMaxSquared) {
            Vector[] forces = potentialCompute.getForces();
            Vector fi = forces[i];
            double fidotdr = fi.dot(dr);
            Vector fj = forces[j];
            double fjdotdr = fj.dot(dr);
            double fdr = fjdotdr - fidotdr;
            for (int k = 0; k < xDataSource.getData().getLength(); k++) {
                double R = xDataSource.getData().getValue(k);
                     double xu = R < r ? 1 : 0; //calcXu(r, u, R);
               //      gSum[k] -= ((xu/(4 * Math.PI)))* wp*beta;               //add once for each atom
                     //replaced fi.dot(dr)-fj.dot(dr) terms with fi.dot(dr)
                    gSum[k] =gSum[k] - ((xu/(4 * Math.PI*r * r * r))-(1/(3*vol)))* fdr*beta/vol;               //add once for each atom
                    //   gSum[k] -= ((-r*r*r/(3*vol)))* wp*beta;               //add once for each atom
                     gSum2[k] = gSum2[k] - ((xu/(4 * Math.PI*r * r * r)) )* fdr*beta/vol;                //add once for each atom
                     thirdterm[k] +=  (beta*fdr/(3*vol*(4*Math.PI*rcforHandfinmap*rcforHandfinmap*rcforHandfinmap/3)));               //add once for each atom

                     newgSum[k]=newgSum[k]-((beta/vol)*-1*fjdotdr*((xu/(4*Math.PI*r * r * r))-(1/(3*vol))));

                     newestgSum[k]=newestgSum[k]+(beta*-1*fjdotdr*(1-xu)/(vol*4*Math.PI*r * r * r));
//
            }
        }
     }

    public double[] getGSum() {
        return gSum;
    }
}
