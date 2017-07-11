/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.statmech;

/**
 * Hard sphere properties based on the Hall EOS.
 */
public class HardSphere {

    protected final double sqrt2 = Math.sqrt(2);

    public HardSphere(){
    }

    /**
     * Returns the compressibility factor for solid (FCC) hard spheres at the
     * given density.
     *
     * @param rho Density of solid hard spheres
     * @return z Compressibility
     */
    public double zSolid(double rho) {
        double beta = 4*(1-rho/sqrt2);
        double z = 12/beta + (2.557696-3) + (0.1253077 + (0.1762393 + (-1.053308 + (2.818621 + (-2.921934 + 1.118413*beta)*beta)*beta)*beta)*beta)*beta;
        return z;
    }

    protected double integrand(double rho){
        double z = zSolid(rho);
        return z/rho;
    }

    /**
     * Returns the difference between absolute free energies between the given
     * densities for solid (FCC) hard spheres.  The result is computed using
     * numeric integration with n points.
     *
     * @param rho1 Density of solid 1
     * @param rho2 Density of solid 2
     * @param n number of points
     * @return Absolute free energy difference
     */
    public double deltaA(double rho1, double rho2, int n){
        double h = (rho2 - rho1)/n;
        double sum = 0.0;
        sum += 0.5*integrand(rho1);
        for (int k=1; k<n; k++) {
            sum += integrand(rho1+k*h);
        }
        sum += 0.5*integrand(rho2);
        return (h*sum);    
    }

    /**
     * Returns the ideal gas free energy for the given density.
     * @param rho Density of ideal gas
     * @return Ideal Gas Free Energy
     */
    public double idFreeEnergy(double rho){return (Math.log(rho)-1.0);}

    public static void main(String[] args) {

        HardSphere hs = new HardSphere();
        double fid = hs.idFreeEnergy(1.04086);
        for(double rho = 1.0; rho < 1.411; rho += 0.01) {
            double deltaF = hs.deltaA(1.04086,rho,100000);
//            System.out.println("rho: "+rho);
//            System.out.println("Ideal Gas Free Energy = "+fid);
//            System.out.println("Excess Free Energy = "+fex);
//            System.out.println("Absolute Free Energy (N= 32) = "+(fex+fid+5.8644));//Fexcess = 5.8644 for N = 32;5.9117 for N = 108 ;5.9208 for N = 256 JCP81 Pg 3191
//            System.out.println("Absolute Free Energy (N=108) = "+(fex+fid+5.9117));//N = 108
//            System.out.println("Absolute Free Energy (N=256) = "+(fex+fid+5.9208));//N = 256
            double bA = (deltaF+fid+5.9208);
            double bmu = hs.deltaA(1.04086,rho,100000)+fid+5.9208 + hs.zSolid(rho);
            System.out.println(String.format("%5.2f  %20.15e  %20.15e", rho, bA, bmu));
        }
    }
}
