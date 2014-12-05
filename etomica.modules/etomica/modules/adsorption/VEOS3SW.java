/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;


public class VEOS3SW extends EOSSW {

    protected double B2, B3;

    protected void reset() {
        double b0 = 2.0/3.0*Math.PI*sigma*sigma*sigma;
        double lambda3 = lambda*lambda*lambda;
        double x = Math.exp(epsilon/temperature)-1;
        B2 = b0*(1 - (lambda3-1)*x);
        B3 = b0*b0*0.125*(5.0 - 17*x - x*x*(48+lambda*lambda-32*lambda3) - x*x*x*(5*lambda3*lambda3-32*lambda3+18*lambda*lambda+26));
    }
    
    public double pressure(double density) {
        System.out.println("z2 z3 "+B2*density+" "+B3*density*density);
        return temperature*density*(1 + B2*density + B3*density*density);
    }
    
    public double dPdRho(double density) {
        return temperature*(1 + 2*B2*density + 3*B3*density*density);
    }
    
    public double mu(double density) {
        System.out.println("mu B2 B3 "+(2*B2*density) +" "+ 1.5*B3*density*density+" "+Math.log(density*temperature));
        return (2*B2*density + 1.5*B3*density*density + Math.log(density))*temperature;
    }
    
    public double rhoForPressure(double pressure) {
        double rho1 = 0;
        double p1 = 0;
        double rho2 = 0.5*pressure/temperature;
        double p2 = pressure(rho2);
        double dpdr = dPdRho(rho2);
        while (p2 < 0 || dpdr < 0) {
            rho2 *= 0.1;
            p2 = pressure(rho2);
            dpdr = dPdRho(rho2);
        }
        boolean first = true;
//        System.out.println("trying for "+pressure);
//        System.out.println("r p "+rho2+" "+p2);
        while (p2 < pressure) {
            double _rho2 = rho2;
            if (dpdr > 0) {
                rho2 += (pressure-p2)/dpdr;
            }
            else {
                rho2 = 0.5*(rho1+rho2);
            }
            rho1 = _rho2;
            p1 = p2;
            if (!first) {
                rho2 *= 1.01;
            }
            first = false;
            double pNew = pressure(rho2);
//            System.out.println("r p "+rho2+" "+pNew);
            if (pNew < p2) {
                throw new RuntimeException("oops");
            }
            p2 = pNew;
            dpdr = dPdRho(rho2);
        }
        double rho = 0;
        while (true) {
            rho = rho1 + (pressure-p1)*(rho2-rho1)/(p2-p1);
//            System.out.println("rho 1 2 new    "+rho1+" "+rho2+" "+rho);
            if (rho == rho1 || rho==rho2) break;
            if (rho < rho1 || rho > rho2) throw new RuntimeException("oops");
            double p = pressure(rho);
            if (p == pressure) break;
            if (p > pressure) {
                p2 = p;
                rho2 = rho;
            }
            else {
                p1 = p;
                rho1 = rho;
            }
        }
        return rho;
    }
    
    public double muForPressure(double pressure) {
        double rho = rhoForPressure(pressure);
        System.out.println("rho p mu "+rho+" "+pressure+" "+mu(rho));
        return mu(rho);
    }
}
