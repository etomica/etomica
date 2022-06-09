/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;
import etomica.units.ElectronVolt;

/**
* LREP-V
 */
public class P2LREPV implements IPotential2 {

    protected double c0, c1, c2, c3, c4, rC;
    protected int m;

    public P2LREPV() {
        this.m = 4;
        this.rC = 6.1;
        this.c0  = ElectronVolt.UNIT.toSim(0.123554);
        this.c1  = ElectronVolt.UNIT.toSim(-0.134361);
        this.c2  = ElectronVolt.UNIT.toSim(0.0543818);
        this.c3  = ElectronVolt.UNIT.toSim(-0.981194E-2);
        this.c4  = ElectronVolt.UNIT.toSim(0.675816E-3);
        System.out.println(" LREP-V parameters:  rC-1: "+ rC +" c0: "+c0+" c1: "+c1+" c2: "+c2+" c3: "+c3+" c4: "+c4);
    }

    /**
     * The energy u.
     */
    public double u(double r2) {
        double r=Math.sqrt(r2);
        if(r<= rC){
            double V1 = (r- rC)*(r- rC)*(r- rC)*(r- rC);
            double V2  = c0 + r*(c1 + r*(c2 + r*(c3 + c4*r))) ;
            return V1*V2;
        }else{
            return 0;
        }
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        double r=Math.sqrt(r2);
        if(r<= rC){
            double V1 = (r- rC)*(r- rC)*(r- rC)*(r- rC);
            double dV1 = m *(r- rC)*(r- rC)*(r- rC);
            double V2  = c0 + r*(c1 + r*(c2 + r*(c3 + c4*r))) ;
            double dV2  = c1+r*(2.0*c2+r*(3.0*c3+4.0*c4*r));
            return r*(V1*dV2+dV1*V2);
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
        if(r<= rC){
            double V1 = (r- rC)*(r- rC)*(r- rC)*(r- rC);
            double dV1 = m *(r- rC)*(r- rC)*(r- rC);
            double ddV1 = m *(m -1.0)*(r- rC)*(r- rC);
            double V2  = c0 + r*(c1 + r*(c2 + r*(c3 + c4*r))) ;
            double dV2  = c1+r*(2.0*c2+r*(3.0*c3+4.0*c4*r));
            double ddV2  = 2.0*c2+r*(6.0*c3+12.0*c4*r);
            return r2*(2.0*dV1*dV2+V1*ddV2+ddV1*V2);
        }else{
            return 0;
        }
    }

    /**
     * Returns the truncation radius.
     */
    public double getRange() { return rC; }


}
