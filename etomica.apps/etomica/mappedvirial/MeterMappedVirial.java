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
import etomica.data.meter.MeterPressure;
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
import etomica.statmech.LennardJones;
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
    protected final double[] c1;
    protected final double[][] cumint;
    protected final MeterPressure meterP;
    protected final AtomPair pair;
    protected final double vol;
    protected final double[] q;
    protected final double[] rCutoff;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final double[] sum, virialSum;
    
    public MeterMappedVirial(ISpace space, IPotentialMaster potentialMaster, Potential2SoftSpherical p2, IBox box, int nbins, double[] rCutoff) {
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
        c1 = new double[rCutoff.length];
        for (int j=0; j<rCutoff.length; j++) {
            c1[j] = Math.log(rCutoff[j]+1)/nbins;
        }
        meterP = new MeterPressure(space);
        meterP.setBox(box);
        meterP.setPotentialMaster(potentialMaster);
        vol = box.getBoundary().volume();
        this.rCutoff = rCutoff;
        q = new double[rCutoff.length];
        sum = new double[rCutoff.length];
        virialSum = new double[rCutoff.length];
        data = new DataDoubleArray(rCutoff.length);
        dataInfo = new DataInfoDoubleArray("mapped virial", Pressure.DIMENSION, new int[]{rCutoff.length});
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    protected double integrand(double r2, int idx) {
        double u = p2.u(r2);
        return Math.exp(-beta*u);
    }

    public void setTemperature(double T) {
        meterP.setTemperature(T);
        beta = 1/T;
        int D = space.D();
        int nbins = cumint[0].length-1;
        for (int j=0; j<rCutoff.length; j++) {
            
            for (int i=1; i<=nbins; i++) {
                double r = Math.exp(c1[j]*i)-1;
                double r2 = r*r;
                double integrand = Math.exp(-beta*p2.u(r2))-1;
                q[j] += (D==2 ? r : r2)*integrand*c1[j]*(r+1);
            }
            q[j] *= 4*Math.PI;

            for (int i=1; i<=nbins; i++) {
                // r = (exp(c1*i) - 1)
                // rMax = (exp(c1*nbins)-1)
                // c1 = ln(xMax+1)/nbins
                double r = Math.exp(c1[j]*i)-1;
                double r2 = r*r;
                cumint[j][i] = cumint[j][i-1] + (D==2 ? r : r2)*integrand(r2, j)*c1[j]*(r+1);
            }
        }

    }
    
    public static void main(String[] args) {
        ISpace space = Space.getInstance(3);
        IBox box = new Box(space);
        double V = 100/0.025;
        double L = Math.pow(V,1.0/3.0);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{L,L,L}));
        P2LennardJones p2 = new P2LennardJones(space);
        P2SoftSphericalTruncated p2t = new P2SoftSphericalTruncated(space, p2, 4.0);
        p2.setBox(box);
        double[] cutoff = new double[]{0.95, 0.97, 1.0, 1.1, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0};
        int nbins = 10000;
        MeterMappedVirial meter = new MeterMappedVirial(space, null, p2t, box, nbins, cutoff);
        meter.setTemperature(1.0);
        for (int i=0; i<cutoff.length; i++) {
            System.out.println(i+" "+cutoff[i]+" "+meter.q[i]);
        }
        int idx = 7;
        double rc = cutoff[idx];
        for (int i=1; i<1000; i++) {
            double r = i*rc/1000;
            double u = p2.u(r*r);
            System.out.println(r+" "+meter.cumint(idx, r, u));
        }
    }

    protected double calcXs(int idx, double r, double u) {
        double qV = q[idx]/vol;
        if (r > rCutoff[idx]) {
            return 3*q[idx]/(4*Math.PI)*(1-qV)/(r*r) - qV*r;
        }
        double y = cumint(idx, r, u);
        return -r + 3*(1-qV)*y*Math.exp(beta*u)/(r*r);
    }

    protected double cumint(int idx, double r, double u) {
        // r = (exp(c1*i) - 1)
        double i = Math.log(r+1)/c1[idx];
        int ii = (int)i;
        double y = cumint[idx][ii] + (cumint[idx][ii+1]-cumint[idx][ii])*(i-ii);
        return y;
    }

    public MyAgent makeAgent(IAtom a) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(MyAgent agent, IAtom atom) {}

    public IData getData() {
        pcForce.reset();
        potentialMaster.calculate(box, allAtoms, pcForce);
        IAtomList list = box.getLeafList();
        for (int j=0; j<sum.length; j++) {
            sum[j] = virialSum[j] = 0;
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
                IVector[] gij = p2.gradient(pair);
                double term = Double.NaN;
                for (int jj=rCutoff.length-1; jj>=0; jj--) {
                    if (r>rCutoff[jj]) {
                        virialSum[jj] += gij[0].dot(dr);
                        continue;
                    }
                    double xs = cumint(jj, r, p2.u(r2));
                    if (Double.isNaN(term)) {
                        term = (fi.dot(dr) + gij[0].dot(dr)) - (fj.dot(dr) + gij[1].dot(dr));
                    }
                    sum[jj] += xs*term;
                }
            }
        }
//        System.out.println("mappedV "+sum[sum.length-1]+" "+B2[B2.length-1]);
        int D = space.D();
        double density = n/vol;
        double[] x = data.getData();

        double vol = box.getBoundary().volume();
        for (int j=0; j<rCutoff.length; j++) {
            x[j] = density/beta - 0.5*q[j]*density*density/beta  - sum[j]/(2*D*vol) + virialSum[j]/(D*vol);
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
