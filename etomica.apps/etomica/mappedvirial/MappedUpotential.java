package etomica.mappedvirial;

 import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomPair;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space3d.Space3D;

/**
 * PotentialCalculation that implements mapped-averaged virial (to get the
 * pressure).
 * 
 * @author Andrew Schultz
 */
public class MappedUpotential implements PotentialCalculation {

    protected final IBox box;
    protected final IteratorDirective allAtoms;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final ISpace space;
    protected double beta;
    protected final IVectorMutable dr;
    protected double c1;
    protected final double[] cumint;
    protected final AtomPair pair;
    protected double vol;
    protected double q;
    protected double qp;
    protected double qp_q;
    protected final int nbins;
    protected double sum;
    protected double x0, vCut;

    public MappedUpotential(ISpace space, IBox box, int nbins, AtomLeafAgentManager<MyAgent> forceManager) {
        this.space = space;
        this.box = box;
        this.nbins = nbins;
        pair = new AtomPair();
        this.forceManager = forceManager;
        allAtoms = new IteratorDirective();
        dr = space.makeVector();
        cumint = new double[nbins+1];
        vol = box.getBoundary().volume();
      //  vol = 1e5;
      }

    public static void main (String[] args) throws IOException{
        Simulation sim = new Simulation(Space3D.getInstance());
        IBox box = new Box(sim.getSpace());
        MappedUpotential pc = new MappedUpotential(sim.getSpace(),box, 1000000, null);
        P2LennardJones potential = new P2LennardJones(sim.getSpace());
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncatedShifted(sim.getSpace(), potential, 4);
        
        pc.setTemperature(2, p2Truncated);
        FileWriter fw = new FileWriter("vb.dat");
        for (int i=13; i<=50; i++) {
            double r = i*0.05;
            if (r>=2.5) r = 2.499999;
         //   double ulrc = potential.uInt(4);
       //     double ulr = potential.uInt(r);
          //  System.out.println(r+" "+pc.calcXs(r, p2Truncated.u(r*r))+" "+(pc.qp/(4*Math.PI*r*r)*(-pc.vol/ pc.q)-((ulrc-ulr)/(r*r)))+(pc.qp_q*r/3));
            fw.write(r+" "+pc.calcXs(r, p2Truncated.u(r*r))+"\n");
        }
        fw.close();
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
            double v = u;
            if(r>vCut) v=0;
            evm1 = Math.exp(-beta*v);
            qp += (D==2 ? r : r2)*evm1*v*c1*(r+1);
            q += (D==2 ? r : r2)*(evm1-1)*c1*(r+1);
            
        }
        qp *= -1*(D==2?2:4)*Math.PI;
        q *= (D==2?2:4)*Math.PI;
        q += vol;
        qp_q = qp/q;

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
            double v = u;
            if(r>vCut) v=0;
            evm1 = Math.exp(-beta*v);
            cumint[i] = cumint[i-1] + (D==2 ? r : r2)*(evm1)*(v+qp_q)*c1*(r+1);
        }
    }
    
    protected double calcXs(double r, double u) {
        double y = cumint(r);
        double v = u;
        if(r>vCut) v=0;
        double evm1 = Math.exp(-beta*v);
        return y/(r*r*evm1);
    }

    protected double cumint( double r) {
        // r = (exp(c1*i) - 1)
        double i = Math.log(r+1)/c1;
        int ii = (int)i;
        double y = cumint[ii] + (cumint[ii+1]-cumint[ii])*(i-ii);
        return y;
    } 

    public double getQP_Q() {
        return qp_q;
    }
    
    public double getqp() {
        return qp;
    }
    
    public void reset() {
        sum = 0;
    }
    
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof Potential2SoftSpherical)) return;
        Potential2SoftSpherical p2 = (Potential2SoftSpherical)potential;
        IAtom a = atoms.getAtom(0);
        IAtom b = atoms.getAtom(1);
//        System.out.println("volume "+vol);
        dr.Ev1Mv2(b.getPosition(),a.getPosition());
        box.getBoundary().nearestImage(dr);
        double r2 = dr.squared();
        double r = Math.sqrt(r2);
        if (r > p2.getRange()) return;
        double fij = p2.du(r2);
        double up = fij/r;
        double vp = up;
        if(r>vCut) sum += -p2.u(r2);
        if (r<x0) {
            IVector fi = forceManager.getAgent(a).force;
            IVector fj = forceManager.getAgent(b).force;
            double u = p2.u(r2);
          //  System.out.println(u+" "+r);
            double fifj = (fi.dot(dr) - fj.dot(dr))/r;
            double xs = calcXs(r, u);    
          //  System.out.println(xs+" "+r);
            double wp = 0.5*fifj;
            sum += xs*beta*(vp-wp);
            
        }
        
    }
    
       public double getPressure() {           
       
        return sum;
        
    }
}
