/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

public abstract class EOSSW {

    protected double sigma, lambda, epsilon, temperature;
    // valid for lambda=1.5, from
    //  "Surface tension and vapor-­liquid phase coexistence of the square-well fluid"
    //    Jayant K. Singh, David A. Kofke, and Jeffrey R. Errington
    //    J. Chem. Phys. (119) 3405
    public static double Tc = 1.2172, rhoc = 0.3079;

    public synchronized double pSat() {
        // lnP vs. 1/T fit for data with lambda=1.5 from 
        // "Surface tension and vapor-­liquid phase coexistence of the square-well fluid"
        //  Jayant K. Singh, David A. Kofke, and Jeffrey R. Errington
        //  J. Chem. Phys. (119) 3405
        // http://link.aip.org/link/doi/10.1063/1.1590313
        if (lambda != 1.5) throw new RuntimeException("oops, data is for lambda=1.5");
        if (temperature > Tc) throw new IllegalArgumentException("supercritical");
        return Math.exp(3.5211 - 7.1939 / temperature);
    }

    public synchronized void setLambda(double newLambda) {
        lambda = newLambda;
        reset();
    }

    public synchronized void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
        reset();
    }

    public synchronized void setSigma(double sigma) {
        this.sigma = sigma;
        reset();
    }

    public synchronized void setTemperature(double newTemperature) {
        temperature = newTemperature;
        reset();
    }

    protected abstract void reset();
    
    public abstract double pressure(double rho);
    
    public abstract double mu(double rho);
    
    public abstract double rhoForPressure(double pressure);
    
    public abstract double muForPressure(double pressure);

}