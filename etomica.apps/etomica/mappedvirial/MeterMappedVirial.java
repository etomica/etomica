package etomica.mappedvirial;

 import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomPair;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.units.Pressure;

public class MeterMappedVirial implements IEtomicaDataSource, AgentSource<IntegratorVelocityVerlet.MyAgent> {

    protected final IPotentialMaster potentialMaster;
    protected final PotentialCalculationForceSum pcForce;
    protected final IBox box;
    protected final IteratorDirective allAtoms;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final ISpace space;
    protected double beta;
    protected final Potential2SoftSpherical p2;
    protected final IVectorMutable dr;
    protected final double c1;
    protected final double[][] cumint;
    protected final AtomPair pair;
    protected final double vol;
    protected final double[] q;
    protected final double[] rCutoff;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final double[] sum;
    protected final double[] epsilon;
    protected final boolean switchev;
    
    public MeterMappedVirial(ISpace space, IPotentialMaster potentialMaster, Potential2SoftSpherical p2, IBox box, int nbins, double[] rCutoff) {
        switchev = false;
        this.space = space;
        this.p2 = p2;
        this.box = box;
        this.potentialMaster = potentialMaster;
        pair = new AtomPair();
        pcForce = new PotentialCalculationForceSum();
        if (box != null) {
            forceManager = new AtomLeafAgentManager<MyAgent>(this, box, MyAgent.class);
            pcForce.setAgentManager(forceManager);
        }
        else {
            forceManager = null;
        }
        allAtoms = new IteratorDirective();
        dr = space.makeVector();
        cumint = new double[rCutoff.length][nbins+1];
        c1 = Math.log(p2.getRange()+1)/nbins;
        vol = box.getBoundary().volume();
        this.rCutoff = rCutoff;
        q = new double[rCutoff.length];
        sum = new double[rCutoff.length];
        data = new DataDoubleArray(rCutoff.length);
        dataInfo = new DataInfoDoubleArray("mapped virial", Pressure.DIMENSION, new int[]{rCutoff.length});
        tag = new DataTag();
        dataInfo.addTag(tag);
        epsilon = new double[rCutoff.length];
        double pc = p2.getRange();
        for (int i=0; i<rCutoff.length; i++) {
            epsilon[i] = (pc-rCutoff[i])*0.5;
            if (epsilon[i] > 0.3*rCutoff[i]) epsilon[i] = 0.3*rCutoff[i];
        }
    }
    
    public double[] getEpsilon() {
        return epsilon;
    }
    
    public void setTemperature(double T) {
        beta = 1/T;
        int D = space.D();
        int nbins = cumint[0].length-1;
        for (int j=0; j<rCutoff.length; j++) {
            
            for (int i=1; i<=nbins; i++) {
                double r = Math.exp(c1*i)-1;
                double r2 = r*r;
                double u = p2.u(r2);
                double evm1 = 0;
                if (switchev) {
                    double eu = Math.exp(-beta*u);
                    evm1 = (eu-1)*(1-theta(r,j));
                }
                else {
                    double v = u*(1-theta(r,j));
                    evm1 = Math.exp(-beta*v)-1;
                }
                q[j] += (D==2 ? r : r2)*evm1*c1*(r+1);
                cumint[j][i] = cumint[j][i-1] + (D==2 ? r : r2)*(evm1+1)*c1*(r+1);
//                if (rCutoff[j]==1 && r>0.864 && r<0.896) System.out.println(r+" "+(evm1+1)+" "+cumint[j][i]);
            }
            q[j] *= 4*Math.PI;
//            System.out.println(rCutoff[j]+" "+q[j]);
        }
    }
    
    public static void main(String[] args) {
        ISpace space = Space.getInstance(3);
        IBox box = new Box(space);
        double V = 100/0.025;
        double L = Math.pow(V,1.0/3.0);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{L,L,L}));
        P2LennardJones p2 = new P2LennardJones(space);
        P2SoftSphericalTruncated p2t = new P2SoftSphericalTruncated(space, p2, 20.0);
        p2.setBox(box);
        double[] cutoff = new double[]{0.95, 0.97, 1.0, 1.1, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0,10.0};
        int nbins = 100000;
        MeterMappedVirial meter = new MeterMappedVirial(space, null, p2t, box, nbins, cutoff);
        meter.setTemperature(1.0);
        for (int i=0; i<cutoff.length; i++) {
            System.out.println(i+" "+cutoff[i]+" "+meter.q[i]);
        }
        double rc = 5; //p2t.getRange();
        for (int idx = 0; idx<cutoff.length; idx++) {
            idx=7;
            for (int i=1; i<1000; i++) {
                double r = i*rc/1000;
                double u = p2.u(r*r);
//                meter.calcXs(idx,r,u);
                System.out.println(r+" "+meter.cumint(idx,r)+" "+meter.calcXs(idx, r, u));
            }
            System.out.println("&");
            break;
        }
    }

    protected double calcXs(int idx, double r, double u) {
        double y = cumint(idx, r);
        double evm1 = 0;
        if (switchev) {
            double eu = Math.exp(-beta*u);
            evm1 = (eu-1)*(1-theta(r,idx));
        }
        else {
            double v = u*(1-theta(r,idx));
            evm1 = Math.exp(-beta*v)-1;
        }
//        System.out.println(r+" "+eu+" "+theta(r,idx)+" "+(evm1+1)+" "+y);
        return -r + space.D()*y/(r*r*(evm1+1));
    }

    protected double cumint(int idx, double r) {
        // r = (exp(c1*i) - 1)
        double i = Math.log(r+1)/c1;
        int ii = (int)i;
        double y = cumint[idx][ii] + (cumint[idx][ii+1]-cumint[idx][ii])*(i-ii);
        return y;
    }

    public MyAgent makeAgent(IAtom a) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(MyAgent agent, IAtom atom) {}
    
    public double theta(double r, int idx) {
        return 0.5*(1+Math.tanh((r-rCutoff[idx])/epsilon[idx]));
    }

    public double dtheta(double r, int idx) {
        double cosh = Math.cosh((r-rCutoff[idx])/epsilon[idx]);
        return 0.5/(epsilon[idx]*cosh*cosh);
    }
    
    public double p(double e, double theta) {
        return 1 + (e-1)*(1-theta);
    }

    public IData getData() {
        pcForce.reset();
        potentialMaster.calculate(box, allAtoms, pcForce);
        IAtomList list = box.getLeafList();
        for (int j=0; j<sum.length; j++) {
            sum[j] = 0;
        }
        int n = list.getAtomCount();
        for (int i=0; i<n; i++) {
            IAtom a = list.getAtom(i);
            IVector fi = forceManager.getAgent(a).force;
            pair.atom0 = a;
            for (int j=i+1; j<n; j++) {
                IAtom b = list.getAtom(j);
                dr.Ev1Mv2(b.getPosition(),a.getPosition());
                box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                double r = Math.sqrt(r2);
                if (r > p2.getRange()) continue;
                IVector fj = forceManager.getAgent(b).force;
                pair.atom1 = b;
                double u = p2.u(r2);
                double eu = Math.exp(-beta*u);
                double fifj = (fi.dot(dr) - fj.dot(dr))/r;
                double fij = p2.du(r2);
                for (int jj=0; jj<rCutoff.length; jj++) {
                    double xs = calcXs(jj, r, u);
                    double theta = theta(r, jj);
                    double up = fij/r;
                    double vp = 0;
                    if (switchev) {
                        double ev = 1 + (eu-1)*(1-theta);
                        vp = (up*eu*(1-theta) + (eu-1)*dtheta(r,jj)/beta)/ev;
                    }
                    else {
                        vp = up*(1-theta) - u*dtheta(r,jj);
                    }
                    double wp = 0.5*fifj;
                    sum[jj] += r*(vp-up) + xs*(vp-wp);
                }
            }
        }
//        System.out.println("mappedV "+sum[sum.length-1]+" "+B2[B2.length-1]);
        int D = space.D();
        double density = n/vol;
        double[] x = data.getData();

        for (int j=0; j<rCutoff.length; j++) {
            x[j] = density/beta - 0.5*q[j]*density*density/beta  + sum[j]/(D*vol);
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
}
