package etomica.normalmode.nptdemo.fluid;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.potential.Potential2Spherical;
import etomica.space.ISpace;

public class Stuff2 implements IStuff {

    protected final IBox box;
    protected final IVectorMutable[] xi;
    protected final IVectorMutable dr;
    protected final int nbins;
    protected final ISpace space;
    protected double rCut;
    protected final double[] cumint;
    protected final double[] sum2;
    protected final double[][] Jij;
    protected final double c1;
    protected final Potential2Spherical p2;
    protected double beta;
    protected boolean doJ = false;
    
    public Stuff2(ISpace space, Potential2Spherical p2, double T, IBox box, double rMax, int nbins) {
        this.space = space;
        this.box = box;
        this.p2 = p2;
        int n = box.getLeafList().getAtomCount();
        xi = new IVectorMutable[n];
        for (int i=0; i<n; i++) {
            xi[i] = space.makeVector();
        }
        Jij = new double[space.D()*n][space.D()*n];
        dr = space.makeVector();
        cumint = new double[nbins];
        sum2 = new double[n];
        c1 = Math.log(rMax+1)/nbins;
        this.nbins = nbins;
        p2.setBox(box);
        setTemperature(T);
    }
    
    public void setDoJacobian(boolean newJ) {
        doJ = newJ;
    }
    
    public void setRCut(double rCut) {
        this.rCut = rCut;
    }

    public void setTemperature(double T) {
        beta = 1/T;
        int D = space.D();
        for (int i=1; i<nbins; i++) {
            // r = (exp(c1*i) - 1)
            // rMax = (exp(c1*nbins)-1)
            // c1 = ln(xMax+1)/nbins
            double r = Math.exp(c1*i)-1;
            double r2 = r*r;
            double u = p2.u(r2);
            cumint[i] = cumint[i-1] + (D==2 ? r : r2)*Math.exp(-beta*u)*c1*(r+1);
            if (cumint[i] < 0) throw new RuntimeException("oops "+i+" "+cumint[i]);
        }
    }
    
    protected double cumint(double r) {
        // r = (exp(c1*i) - 1)
        double i = Math.log(r+1)/c1;
        int ii = (int)i;
        double y = cumint[ii] + (cumint[ii+1]-cumint[ii])*(i-ii);
        if (y<0) throw new RuntimeException("oops "+i+" "+y+" "+ii+" "+cumint[ii]+" "+cumint[ii+1]);
        return y;
    }

    
    public IVector[] stuff() {
        IAtomList atoms = box.getLeafList();
        int n = atoms.getAtomCount();
        for (int i=0; i<n; i++) {
            xi[i].E(0);
            sum2[i] = 0;
        }
        IBoundary boundary = box.getBoundary();
        double vol = boundary.volume();
        int D = space.D();
        double volFac = 1.0/(vol*D);
        double m = 2;
        double cp = 3;
        for (int i=0; i<n; i++) {
            IAtom atomi = atoms.getAtom(i);
            IVector ri = atomi.getPosition();
            for (int j=i+1; j<n; j++) {
                dr.Ev1Mv2(ri, atoms.getAtom(j).getPosition());
                boundary.nearestImage(dr);
                double r2 = dr.squared();
                if (r2 > rCut*rCut) continue;

                double r = Math.sqrt(r2);
                double cumintr = cumint(r);
                double fac = 0.5*volFac/(1+Math.pow(cumintr,cp));
//                System.out.println(r+" "+volFac+" "+0.5/(1+Math.pow(cumint(r),3)));
//                System.out.println(r+" "+cumint(r));
                double fac2 = Math.pow(cumintr, -m);
//                if (i==0) System.out.println(j+" "+dr+" "+r+" "+fac);
//                System.out.println(r+" "+r/volFac+" "+r*fac);
                xi[i].PEa1Tv1(-fac*fac2, dr);
                xi[j].PEa1Tv1(fac*fac2, dr);
                
                sum2[i] += fac2;
                sum2[j] += fac2;
                
                if (doJ) {
                    double dSum11 = fac*fac2;
                    double e = Math.exp(-beta*p2.u(r2));
                    double dSum12 = fac*fac*fac2*fac2 * (2/volFac*Math.pow(cumintr,cp-1)*(Math.pow(r,D-1)*e*(r-1)+cumintr*m)/r*Math.pow(r-1,m-1));
                    for (int k=0; k<D; k++) {
                        for (int l=0; k<D; l++) {
                            Jij[D*i+k][D*j+l] = dSum11 - dr.getX(k)*dr.getX(l)*dSum12;
                            Jij[D*i+l][D*j+k] = Jij[D*i+k][D*j+l];
                            Jij[D*j+k][D*i+l] = Jij[D*i+k][D*j+l];
                            Jij[D*j+l][D*i+k] = Jij[D*i+k][D*j+l];
                        }
                    }
                }
            }
//            if (i==0) System.out.println(xi2[i]);
        }
        for (int i=0; i<n; i++) {
            IAtom atomi = atoms.getAtom(i);
            IVector ri = atomi.getPosition();

            if (sum2[i] > 0) {
                xi[i].TE(1/sum2[i]);
            }

            if (doJ) {
                for (int k=0; k<D; k++) {
                    for (int l=0; l<D; l++) {
                        Jij[D*i+k][D*i+l] = 0;
                        for (int j=0; j<n; j++) {
                            if (j==i) continue;
                            Jij[D*i+k][D*i+l] += Jij[D*i+k][D*j+l];
                        }
                    }
                }
            }
//            xi[i].PEa1Tv1(volFac, ri);
        }
        return xi;
    }
    
    public double getJ() {
        //XXX this method doesn't work.  needs eigenvalues
        if (!doJ) throw new RuntimeException("You need to call setDoJacobian!");
        org.apache.commons.math3.linear.Array2DRowRealMatrix J = new Array2DRowRealMatrix(Jij);
        org.apache.commons.math3.linear.LUDecomposition lud = new org.apache.commons.math3.linear.LUDecomposition(J);
        return lud.getDeterminant();
    }
}
