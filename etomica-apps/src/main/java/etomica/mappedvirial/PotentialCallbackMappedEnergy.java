/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;

import etomica.atom.AtomPair;
import etomica.box.Box;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2Soft;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;

import java.io.IOException;

/**
 * PotentialCalculation that implements mapped-averaged framework (to get the energy).
 *
 * @author Akshara Goyal
 */
public class PotentialCallbackMappedEnergy implements PotentialCallback {

    protected final PotentialCompute potentialMaster;
    protected final Box box;
    protected double beta;
    protected final Vector dr;
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
    protected double vShift;

    public PotentialCallbackMappedEnergy(PotentialCompute potentialMaster, Box box, int nbins) {
        this.potentialMaster = potentialMaster;
        this.box = box;
        this.nbins = nbins;
        pair = new AtomPair();
        dr = box.getSpace().makeVector();
        cumint = new double[nbins+1];
        vol = box.getBoundary().volume();
    }

    public static void main (String[] args) throws IOException{
        Simulation sim = new Simulation(Space3D.getInstance());
        Box box = new Box(sim.getSpace());

        PotentialCallbackMappedEnergy pc = new PotentialCallbackMappedEnergy(null, box, 1000000);
        P2LennardJones potential = new P2LennardJones();
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(potential, 4);
        //  System.out.println(vol1);
        pc.setVolume(99999.99999999997);
        pc.setTemperature(2, p2Truncated);
       // FileWriter fw = new FileWriter("vb.dat");
        for (int i=10; i<45; i++) {
            double r = i*0.1;
            if (r>=4) r = 3.99999999;
            System.out.println(r+" "+pc.calcXs(r, p2Truncated.u(r*r))+" ");
         //   fw.write(r+" "+pc.calcXs(r, p2Truncated.u(r*r))+"\n");
        }

       // fw.close();
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

    public void setTemperature(double T, Potential2Soft p2) {
        beta = 1/T;
        double rc = p2.getRange();
        x0 = rc;
        q = 0;
        qp = 0;
        qp_q = 0;

        if (vCut==0) vCut = x0;
        vShift = -p2.u(vCut*vCut)+0.0494908;
        c1 = Math.log(rc+1)/nbins;
        int D = box.getSpace().D();
        for (int i=1; i<=nbins; i++) {
            double r = Math.exp(c1*i)-1;
            double r2 = r*r;
            if (r >= p2.getRange() ) {
                if (i<nbins) throw new RuntimeException("oops "+i+" "+nbins+" "+r+" "+p2.getRange());
                r = Math.exp(c1*i*0.9999)-1;
                r2 = r*r;
            }
            double u = p2.u(r2);
            double v = calcV(r,u);
            double evm1 = Math.exp(-beta*v);
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
            double v = calcV(r,u);
            double evm1 = Math.exp(-beta*v);
            cumint[i] = cumint[i-1] + (D==2 ? r : r2)*(evm1)*(v+qp_q)*c1*(r+1);
        }

    }

    protected double calcXs(double r, double u) {
        double y = cumint(r);
        double v = calcV(r,u);
        double evm1 = Math.exp(-beta*v);
        return y/(r*r*evm1);
    }

    protected double calcV(double r,double u){

        if(r>vCut)
            return 0;
        return u+vShift;
    }

    protected double cumint( double r) {
        // r = (exp(c1*i) - 1)
        double i = Math.log(r+1)/c1;
        int ii = (int)i;
        return cumint[ii] + (cumint[ii+1]-cumint[ii])*(i-ii);
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

    @Override
    public void pairCompute(int i, int j, Vector dr, double[] u012) {
        double r2 = dr.squared();
        double r = Math.sqrt(r2);
        double fij = u012[1];
        double vp = fij/r;
        double u= u012[0];
        double v = calcV(r,u);
        sum += v-u;

        if (r<x0) {
            Vector[] forces = potentialMaster.getForces();
            Vector fi = forces[i];
            Vector fj = forces[j];
            double fifj = (fi.dot(dr) - fj.dot(dr))/r;
            double xs = calcXs(r, u);
            double wp = 0.5*fifj;
            sum += xs*beta*(vp-wp);
        }
    }

    public double getEnergy() {

        return sum;

    }

}
