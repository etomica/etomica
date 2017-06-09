package etomica.mappedvirial;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;


/**
 * PotentialCalculation that implements mapped-averaged virial (to get the
 * pressure).
 * 
 * @author Andrew Schultz
 */
public class PotentialCalculationMappedVirialV implements PotentialCalculation {

    protected final Box box;
    protected final IteratorDirective allAtoms;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final Space space;
    protected double beta;
    protected final Vector dr;
    protected double c1;
    protected final double[] cumint;
    protected final AtomPair pair;
    protected double vol;
    protected double q;
    protected double qu;
    protected double vShift;
    protected final int nbins;
    protected double sum;
    protected double x0, vCut;
    protected VFunc vf;
    protected boolean customVF;

    public PotentialCalculationMappedVirialV(Space space, Box box, int nbins, AtomLeafAgentManager<MyAgent> forceManager) {
        this.space = space;
        this.box = box;
        this.nbins = nbins;
        pair = new AtomPair();
        this.forceManager = forceManager;
        allAtoms = new IteratorDirective();
        dr = space.makeVector();
        cumint = new double[nbins+1];
        vol = box.getBoundary().volume();
    }

    public static void main2(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
        Box box = new Box(sim.getSpace());
        PotentialCalculationMappedVirialV pc = new PotentialCalculationMappedVirialV(sim.getSpace(),box, 1000000, null);
        P2LennardJones potential = new P2LennardJones(sim.getSpace());
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncatedForceShifted(sim.getSpace(), potential, 2.5);
        pc.setTemperature(1, p2Truncated);
        for (int i=10; i<30; i++) {
            double r = i*0.1;
            if (r>=2.5) r = 2.499999;
            System.out.println(r+" "+pc.calcXs(r));
        }
    }
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(space, potential, 4);
        USine uSine = new USine(p2, new double[][]{{0.87624,2.17775},{0.684801,0.647015},{-0.748087,1.1296},{0.692743,-1.07915}}, 2.0);
        for (int i=80; i<400; i++) {
            double r = i*0.01;
            System.out.println(r+" "+uSine.v(r*r,r)+" "+uSine.vp(r*r,r)+" "+p2.u(r*r)+" "+p2.du(r*r)/r);
        }
    }

    public double getX0() {
        return x0;
    }
    
    public void setVCut(double newVCut) {
        vCut = newVCut;
    }
    
    public void setVFunc(VFunc newVF) {
        vf = newVF;
    }

    /**
     * Sets volume to an arbitrary value (instead of the box volume).  This is
     * used for 1/(1+q/V) terms.
     */
    public void setVolume(double newVol) {
        vol = newVol;
    }

    public void setTemperature(double T, P2SoftSphericalTruncated p2) {
        beta = 1/T;
        double rc = p2.getRange();
        x0 = rc*0.95;
        if (vCut==0) vCut = x0;
        c1 = Math.log(rc+1)/nbins;
        int D = space.D();
        qu = 0;
        if (vf==null) {
            UShift ushift = new UShift(p2);
            ushift.setVCut(vCut);
            vf = ushift;
        }
        for (int i=1; i<=nbins; i++) {
            double r = Math.exp(c1*i)-1;
            double r2 = r*r;
            if (r >= p2.getRange()) {
                if (i<nbins) throw new RuntimeException("oops "+i+" "+nbins+" "+r+" "+p2.getRange());
                r = Math.exp(c1*i*0.9999)-1;
                r2 = r*r;
            }
            double u = p2.u(r2);
            double evm1 = 0;
            double v = vf.v(r2,r);
            evm1 = Math.exp(-beta*v)-1;
            q += (D==2 ? r : r2)*evm1*c1*(r+1);
            double eum1 = Math.exp(-beta*u)-1;
            qu += (D==2 ? r : r2)*eum1*c1*(r+1);
            cumint[i] = cumint[i-1] + (D==2 ? r : r2)*(evm1+1)*c1*(r+1);
        }
        q *= (D==2?2:4)*Math.PI;
        qu *= (D==2?2:4)*Math.PI;
    }
    
    protected double calcXs(double r) {
        double y = cumint(r);
        double v = vf.v(r*r,r);
        double evm1 = Math.exp(-beta*v)-1;
        return -r + space.D()/(1+q/vol)*y/(r*r*(evm1+1));
    }

    protected double cumint( double r) {
        // r = (exp(c1*i) - 1)
        double i = Math.log(r+1)/c1;
        int ii = (int)i;
        double y = cumint[ii] + (cumint[ii+1]-cumint[ii])*(i-ii);
        return y;
    }

    public double getQ() {
        return q;
    }
    
    public double getQU() {
        return qu;
    }

    public void reset() {
        sum = 0;
    }
    
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof Potential2SoftSpherical)) return;
        Potential2SoftSpherical p2 = (Potential2SoftSpherical)potential;
        IAtom a = atoms.getAtom(0);
        IAtom b = atoms.getAtom(1);

        dr.Ev1Mv2(b.getPosition(),a.getPosition());
        box.getBoundary().nearestImage(dr);
        double r2 = dr.squared();
        double r = Math.sqrt(r2);
        if (r > p2.getRange()) return;
        double fij = p2.du(r2);
        double up = fij/r;
        double vp = vf.vp(r2,r);
        sum += r*(vp-up);
        if (r<x0) {
            Vector fi = forceManager.getAgent(a).force;
            Vector fj = forceManager.getAgent(b).force;
            double fifj = (fi.dot(dr) - fj.dot(dr))/r;
            double xs = calcXs(r);
            double wp = 0.5*fifj;
            sum += xs*(vp-wp);
//            System.out.println(r+" "+xs+" "+up+" "+vp);
        }
    }

    public double getPressure() {
        int D = space.D();
        return sum/(D*box.getBoundary().volume());
    }
    
    public interface VFunc {
        public double v(double r2, double r);
        public double vp(double r2, double r);
    }
    
    public static class UShift implements VFunc {

        protected final P2SoftSphericalTruncated p2;
        protected double vShift, vCut2;
        
        public UShift(P2SoftSphericalTruncated p2) {
            this.p2 =p2;
            vCut2 = 0.95*p2.getRange();
            vShift = 0;
            vShift = -v(vCut2, Math.sqrt(vCut2));
        }
        
        public void setVCut(double newVCut) {
            vCut2 = newVCut*newVCut;
            vShift = -p2.u(vCut2);
        }
        
        public double v(double r2, double r) {
            if (r2>vCut2) return 0;
            return p2.u(r2)+vShift;
        }

        public double vp(double r2, double r) {
            if (r2>vCut2) return 0;
            return p2.du(r2)/r;
        }
    }
    
    public static class USine implements VFunc {
        protected final P2SoftSphericalTruncated p2;
        protected double vShift, vCut2;
        protected double[][] c;
        protected double beta;
        protected double epsilon = 0.02;

        public USine(P2SoftSphericalTruncated p2, double[][] sineParameters, double temperature) {
            this.p2 = p2;
            c = sineParameters;
            this.beta = 1/temperature;
            vCut2 = 0.95*p2.getRange();
            vCut2 *= vCut2;
            vShift = 0;
            vShift = -v(vCut2, Math.sqrt(vCut2));
        }
        
        public void setVCut(double newVCut) {
            vCut2 = newVCut*newVCut;
            vShift = -v(vCut2,Math.sqrt(vCut2));
        }
        
        public double v(double r2, double r) {
            if (r2>vCut2) return 0;
            double f = c[0][0]-1.01;
            double a = c[0][1];
            for (int i=1; i<c.length; i++) {
                f += c[i][0]*Math.sin(i*a*r+c[i][1]);
            }
//            System.out.print(r+" "+f);
            f *= 0.5*(1+Math.tanh((r-1)/epsilon));
            f += 1;
//            System.out.println(" "+f);
            double e = Math.exp(-beta*p2.u(r2));
            return -Math.log(1+(e-1)*f)/beta + vShift;
        }

        public double vp(double r2, double r) {
            if (r2>vCut2) return 0;
            // v = -ln(1+(e-1)*f)/beta + vShift
            // dv = -(f*de + (e-1)*df)/(beta*x)   x=1+(1-e)*f
            // de = e*(-du/T)
            // df = sum [ci0 cos(a*r + ci1)] 
            double f0 = c[0][0]-1;
            double df0 = 0;
            double a = c[0][1];
            for (int i=1; i<c.length; i++) {
                double sin = c[i][0]*Math.sin(i*a*r+c[i][1]);
                f0 += sin;
                df0 += c[i][0]*Math.cos(i*a*r+c[i][1]) + i*a*sin;
            }
            double theta = 0.5*(1+Math.tanh((r-1)/epsilon));
            double f = f0 * theta;
            f += 1;
            double e = Math.exp(-beta*p2.u(r2));
            double de = -e*p2.du(r2)/r*beta;
            double cosh = Math.cosh((r-1)/epsilon);
            double dTheta = 0.5/(0.1*cosh*cosh);
            double df = f0*dTheta + df0*theta;
            
            return -(f*de + (e-1)*df)/(1+(e-1)*f)/beta;
        }
    }
}
