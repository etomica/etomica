/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

/**
 * EOS for square-well fluids as described by
 * "Generalized equation of state for square-well potentials of variable range"
 *   B. H. Patel, H. Docherty, S. Varga, A. Galindo and G. C. Maitland
 *   Mol. Phys. 2005 (103) 129
 * 
 * @author Andrew Schultz
 */
public class EOSSWPatel extends EOSSW {

    protected double c1, c2, c3;
    protected double alpha;
    
    protected void reset() {
        // we care about lambda=1.5, so use the simple formulation from the older paper.
        c1 = 2.25855-1.50349*lambda+0.249434*lambda*lambda;
        c2 = -0.669270+1.40049*lambda-0.827739*lambda*lambda;
        c3 = 10.1576-15.0427*lambda+5.30827*lambda*lambda;
        alpha = 4*(lambda*lambda*lambda-1);
    }
    
    public double A1(double eta) {
        double etaEff = c1*eta+c2*eta*eta+c3*eta*eta*eta;
        double onemetaEff = 1-etaEff;
        return -epsilon*eta*alpha*(1-etaEff/2)/(onemetaEff*onemetaEff*onemetaEff);
    }
    
    public double A(double density) {
        double eta = density*Math.PI/6*sigma*sigma*sigma;
        double AHS = temperature*(4*eta-3*eta*eta)/((1-eta)*(1-eta));
        double A1 = A1(eta);
        double onemeta = 1-eta;
        double KHS = onemeta*onemeta*onemeta*onemeta/(1+4*eta*(1+eta));
        double A2 = epsilon*0.5*KHS*eta*(A1(eta*1.0001)-A1(eta*0.9999))/(eta*0.0002)/temperature;
        return temperature*(Math.log(density)-1) + AHS+A1+A2;
    }
    
    public double pressure(double density) {
        double dAdrho = (A(density*1.0001)-A(density*0.9999))/(density*0.0002);
        return density*density*dAdrho;
    }

    public double dPdRho(double density) {
        return (pressure(density*1.0001)-pressure(density*0.9999))/(density*0.0002);
    }

    public double mu(double density) {
        double A = A(density);
        double P = pressure(density);
        return A + P/density;
    }

    public double rhoForPressure(double pressure) {
//        if (temperature>0) {
//            System.out.println("start trying for "+pressure+" at T="+temperature);
//        }
        double rho1 = 0;
        double p1 = 0;
        double rho2 = 0.5*pressure/temperature;
        double p2 = pressure(rho2);
        double dpdr = dPdRho(rho2);
        if (temperature < Tc && pressure > pSat()) {
            //liquid
            rho2 = rhoc;
            dpdr = dPdRho(rho2);
            while (dpdr < 0 || p2 < 0) {
                rho2 *= 1.1;
                p2 = pressure(rho2);
                dpdr = dPdRho(rho2);
            }
            p2 = pressure(rho2);
            if (p2 < pressure) {
                // p2 is too low, use it for p1, and we'll find p2 below
                rho1 = rho2;
                p1 = p2;
                rho2 *= 1.1;
                p2 = pressure(rho2);
                dpdr = dPdRho(rho2);
            }
            else {
                // p2 is too high, which is OK, but we need to find a rho1<rho
                // step backwards using dpdr until we find a p1.
                rho1 = rho2;
                p1 = p2;
                boolean first = true;
                while (p1 > pressure) {
                    double _rho1 = rho1;
                    if (dpdr > 0) {
                        rho1 += (pressure-p1)/dpdr;
                    }
                    else {
                        // overstepped (backwards)
                        rho1 = 0.5*(rho1+rho2);
                    }
                    rho2 = _rho1;
                    p2 = p1;
                    if (!first) {
                        rho1 *= 0.999;
                    }
                    first = false;
                    double pNew = pressure(rho1);
//                    System.out.println("r p "+rho2+" "+pNew);
                    if (pNew > p2) {
                        throw new RuntimeException("oops");
                    }
                    p1 = pNew;
                    dpdr = dPdRho(rho1);
                }
            }
        }
        while (p2 < 0 || dpdr < 0) {
            rho2 *= 0.1;
            p2 = pressure(rho2);
            dpdr = dPdRho(rho2);
        }
        boolean first = true;
//        if (temperature>0) {
//            System.out.println("trying for "+pressure+" at T="+temperature);
//            System.out.println("r p "+rho2+" "+p2);
//        }
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
//            if (temperature>0) System.out.println("r p "+rho2+" "+pNew+" "+temperature);
            if (pNew < p2) {
                throw new RuntimeException("oops");
            }
            p2 = pNew;
            dpdr = dPdRho(rho2);
        }
        double rho = 0;
//        System.out.println(rho1+" "+p1);
//        System.out.println(rho2+" "+p2);
        while (true) {
            rho = rho1 + (pressure-p1)*(rho2-rho1)/(p2-p1);
//            System.out.println("rho 1 2 new    "+rho1+" "+rho2+" "+rho);
            if (rho == rho1 || rho==rho2) break;
            if (rho < rho1 || rho > rho2) throw new RuntimeException("oops");
            double p = pressure(rho);
//            System.out.println(rho+" "+p);
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
//        System.out.println("rho p mu "+rho+" "+pressure+" "+mu(rho));
        return mu(rho);
    }
}
