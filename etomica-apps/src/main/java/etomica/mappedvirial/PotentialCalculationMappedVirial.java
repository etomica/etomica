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

//import etomica.potential.P2SoftSphericalTruncatedForceShifted;


/**
 * PotentialCalculation that implements mapped-averaged virial (to get the
 * pressure).
 * 
 * @author Andrew Schultz
 */
public class PotentialCalculationMappedVirial implements PotentialCalculation {

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

    public PotentialCalculationMappedVirial(Space space, Box box, int nbins, AtomLeafAgentManager<MyAgent> forceManager) {
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

    public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
        Box box = new Box(sim.getSpace());
        PotentialCalculationMappedVirial pc = new PotentialCalculationMappedVirial(sim.getSpace(),box, 1000000, null);
        P2LennardJones potential = new P2LennardJones(sim.getSpace());
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(sim.getSpace(), potential, 2.5);
        pc.setTemperature(1, p2Truncated);
        for (int i=10; i<30; i++) {
            double r = i*0.1;
            if (r>=2.5) r = 2.499999;
            System.out.println(r+" "+pc.calcXs(r, p2Truncated.u(r*r)));
        }
    }
    
    public double getX0() {
        return x0;
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

    public void setTemperature(double T, Potential2SoftSpherical p2) {
        beta = 1/T;
        double rc = p2.getRange();
        x0 = rc*0.95;
        if (vCut==0) vCut = x0;
        c1 = Math.log(rc+1)/nbins;
        int D = space.D();
        qu = 0;
        vShift = -p2.u(vCut*vCut);
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
            double v = u + vShift;
            if (r>vCut) v = 0;
            evm1 = Math.exp(-beta*v)-1;
            q += (D==2 ? r : r2)*evm1*c1*(r+1);
            double eum1 = Math.exp(-beta*u)-1;
            qu += (D==2 ? r : r2)*eum1*c1*(r+1);
            cumint[i] = cumint[i-1] + (D==2 ? r : r2)*(evm1+1)*c1*(r+1);
        }
        q *= (D==2?2:4)*Math.PI;
        qu *= (D==2?2:4)*Math.PI;
    }
    
    protected double calcXs(double r, double u) {
        double y = cumint(r);
        double v = u + vShift;
        if (r > vCut) v = 0;
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
        double vp = up;
        if (r>vCut) vp = 0;
        sum += r*(vp-up);
        if (r<x0) {
            Vector fi = forceManager.getAgent(a).force;
            Vector fj = forceManager.getAgent(b).force;
            double u = p2.u(r2);
            double fifj = (fi.dot(dr) - fj.dot(dr))/r;
            double xs = calcXs(r, u);
            double wp = 0.5*fifj;
            sum += xs*(vp-wp);
        }
    }

    public double getPressure() {
        int D = space.D();
        return sum/(D*box.getBoundary().volume());
    }
}
