/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.ElectronVolt;

/**
* LREP-Phi
 */
public class P2LREPPhi implements IPotential2 {

    protected double rC, r0, A, B;
    protected int n;

    public P2LREPPhi() {
        this.n = 6;
        this.rC = 7.8;
        this.A   = ElectronVolt.UNIT.toSim(Math.sqrt(0.656618E-4));
        this.A = A*A;
        this.B   = 1.836569;
        this.r0  = 2.552655;
        System.out.println(" LREP-Phi parameters: rC-2: "+ rC + " A: " + A + " B: " + B + " r0: " + r0);
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        double r=Math.sqrt(r2);
        if(r<=rC){
            double Phi1 = (r-rC)*(r-rC)*(r-rC)*(r-rC)*(r-rC)*(r-rC);
            double Phi2 = Math.exp(-B*(r/r0-1.0));
            return A*Phi1*Phi2;
        }else{
            return 0;
        }
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double r=Math.sqrt(r2);
        if(r<=rC){
            double Phi1 = (r-rC)*(r-rC)*(r-rC)*(r-rC)*(r-rC)*(r-rC);
            double dPhi1 = n *(r-rC)*(r-rC)*(r-rC)*(r-rC)*(r-rC);
            double Phi2 = Math.exp(-B*(r/r0-1.0));
            double dPhi2 = -B/r0*Phi2;
            return r*A*(Phi1*dPhi2+dPhi1*Phi2);
        }else{
            return 0;
        }
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        double r=Math.sqrt(r2);
        if(r<=rC){
            double Phi1 = (r-rC)*(r-rC)*(r-rC)*(r-rC)*(r-rC)*(r-rC);
            double dPhi1 = n *(r-rC)*(r-rC)*(r-rC)*(r-rC)*(r-rC);
            double ddPhi1 = n *(n -1.0)*(r-rC)*(r-rC)*(r-rC)*(r-rC);
            double Phi2 = Math.exp(-B*(r/r0-1.0));
            double dPhi2 = -B/r0*Phi2;
            double ddPhi2 = (B/r0)*(B/r0)*Phi2;
            return r2*A*(2.0*dPhi1*dPhi2+Phi1*ddPhi2+ddPhi1*Phi2);
        }else{
            return 0;
        }
    }

    /**
     * Returns the truncation radius.
     */
    public double getRange() { return rC; }

}
