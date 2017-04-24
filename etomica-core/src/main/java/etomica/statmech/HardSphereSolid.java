/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.statmech;

import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Hard sphere solid EOS.  Free energy at one density taken from
 * http://dx.doi.org/10.1063/1.4818990
 * 
 * Z fit to simulation data for perfect FCC crystal.
 * 
 * @author Andrew Schultz
 */
public class HardSphereSolid {

    protected static final double sqrt2 = Math.sqrt(2);
    // fitting parameters for Z vs. x
    protected static final double[] c = new double[]{-4.338711486444215e-01,3.908120499110030e-01,5.954516832360230e-01,1.859059648880398e+01,-2.881675745716088e+02,3.206943027086472e+03,-2.103760905875148e+04,8.317792071175583e+04,-1.806368017944632e+05,1.702057910137580e+05};
    protected static final double rho0 = 1.0409;
    protected static final double A0 = 5.91901 + idFreeEnergy(rho0);  // Helmholtz free energy at rho=1.0409, A=5.91901(3)
    protected static final double[] cdz = new double[]{4.664021897071290e+00,3.583111436885733e+02,-7.251681562320791e+03,2.057933864384254e+05,1.356484539042134e+06};
    protected static final double[] cc = new double[]{6.125493227538938e-01,2.458082847623682e-01,9.725067704701637e+00,-6.145649875550352e+01,1.631732416910175e+02};

    public HardSphereSolid(){
    }

    public static double calcdMudP(double rho, double drho) {
        double rhop = rho + drho;
        double rhom = rho - drho;
        double Zp = zSolid(rhop);
        double Pp = Zp*rhop;
        double Ap = Asolid(rhop, 100000);
        double mup = Ap + Zp;
        double Zm = zSolid(rhom);
        double Pm = Zm*rhom;
        double Am = Asolid(rhom, 100000);
        double mum = Am + Zm;
        return (mup - mum)/(Pp - Pm);
    }

    public static double calcdPdrho(double rho, double drho) {
        double rhop = rho + drho;
        double rhom = rho - drho;
        double Zp = zSolid(rhop);
        double Pp = Zp*rhop;
        double Zm = zSolid(rhom);
        double Pm = Zm*rhom;
        return (Pp - Pm)/(rhop - rhom);
    }

    
    public static double calcdmudrho(double rho, double drho) {
        double rhop = rho + drho;
        double rhom = rho - drho;
        double Zp = zSolid(rhop);
        double Ap = Asolid(rhop);
        double mup = Zp + Ap;
        double Zm = zSolid(rhom);
        double Am = Asolid(rhom);
        double mum = Zm + Am;
        return (mup - mum)/(rhop - rhom);
    }

    /**
     * Returns the compressibility factor for solid (FCC) hard spheres at the
     * given density.
     * @param rho Density
     * @return Compressibility
     *
     */
    public static double zSolid(double rho) {
        double x = (1-rho/sqrt2);
        return 3/x + zSolidCorrection(rho);
    }

    /**
     * Returns the correction to high density expression for the
     * compressibility factor for solid (FCC) hard spheres at the given
     * density.
     * @param rho Density
     * @return high density correction to the compressibility
     */
    public static double zSolidCorrection(double rho) {
        double x = (1-rho/sqrt2);
        double z = 0;
        double xp = 1;
        for  (int i=0; i<c.length; i++) {
            z += c[i]*xp;
            xp *= x;
        }
        return z;
    }
    
    public static double densityForPressure(double pressure) {
        if (pressure < 10.24) throw new RuntimeException("looks like your pressure corresponds to density<1");
        double rho1 = 2*pressure/(6+sqrt2*pressure);
        double p1 = zSolid(rho1)*rho1;
        if (p1 == pressure) return rho1;
        double p2 = pressure-(p1-pressure)*3;
        double rho2 = 2*p2/(6+sqrt2*p2);
        if (rho2 == rho1) throw new RuntimeException("oops");
        p2 = zSolid(rho2)*rho2;
        if (p2 == pressure) return rho2;
        if ((p1-pressure)*(p2-pressure) > 0) {
            System.out.println("target: "+pressure);
            System.out.println(rho1+" "+p1);
            System.out.println(rho2+" "+p2);
            throw new RuntimeException("oops");
        }
        if (rho2 < rho1) {
            double t = rho1;
            rho1 = rho2;
            rho2 = t;
            t = p1;
            p1 = p2;
            p2 = t;
        }
        while (true) {
            double rho3 = (rho1+rho2)/2;
            if (rho3 == rho1) {
                return rho1;
            }
            if (rho3 == rho2) {
                return rho2;
            }
            double p3 = zSolid(rho3)*rho3;
            if (p3 == pressure) return rho3;
            if (p3 == p1 || p3 == p2) return rho3;
            if (p3 < p1 || p3 > p2) throw new RuntimeException("oops");
            if (p3 > pressure) {
                rho2 = rho3;
                p2 = p3;
            }
            else {
                rho1 = rho3;
                p1 = p3;
            }
        }
    }

    public static double densityForMu(double mu) {
        if (mu < 14.77) throw new RuntimeException("looks like your mu corresponds to density<1");
        double rho1 = 1.00;
        double rho2 = Math.sqrt(2);
        double mu1 = zSolid(rho1) + Asolid(rho1, 100000);
        double mu2 = Double.POSITIVE_INFINITY;
        while (true) {
            double rho3 = (rho1+rho2)/2;
            if (rho3 == rho1) {
                return rho1;
            }
            if (rho3 == rho2) {
                return rho2;
            }
            double mu3 = zSolid(rho3)*rho3;
            if (mu3 == mu) return rho3;
            if (mu3 == mu1 || mu3 == mu2) return rho3;
            if (mu3 < mu1 || mu3 > mu2) throw new RuntimeException("oops");
            if (mu3 > mu) {
                rho1 = rho3;
                mu1 = mu3;
            }
            else {
                rho2 = rho3;
                mu2 = mu3;
            }
        }
    }

    /**
     * Returns the difference between absolute free energies between the given
     * densities for solid (FCC) hard spheres.  The result is computed using
     * numeric integration with n points.
     * @param rho Density
     * @return Absolute free energy difference
     */
    public static double Asolid(double rho){
        return Asolid(rho, 10000);
    }
    
    public static double Asolid(double rho, int n) {
        // integral from rho0 to rho of 3/((1-rho/sqrt2)*rho)
        // this is the high-density expression for the change in the free energy
        double dA0 = -3*(Math.log((1-rho/sqrt2)/(1-rho0/sqrt2)*rho0/rho));
        // now use trapezoid rule to integrate the correction
        double h = (rho - rho0)/n;
        double sum = 0.0;
        sum += 0.5*zSolidCorrection(rho0)/rho0;
        for (int k=1; k<n; k++) {
            sum += zSolidCorrection(rho0+k*h)/(rho0+k*h);
        }
        sum += 0.5*zSolidCorrection(rho)/rho;
        return A0 + dA0 + h*sum;
    }

    /**
     * Returns the fraction of lattice sites with a vacancy for the given
     * lattice density.
     * @param rhoLat lattice density
     * @return Vacant fraction of lattice sites
     */
    public double getVacancyFraction(double rhoLat) {
        double x = 1 - rhoLat/Math.sqrt(2);
        double fitSum = 0;
        double xp = 1;
        for (int i=0; i<cc.length; i++) {
          fitSum += cc[i]*xp;
          xp *= x;
        }
        // fitSum = (ln(c)/Zlat + 1)/x
        double c = (fitSum*x-1)*zSolid(rhoLat);
        return c;
    }

    /**
     * Returns the difference between the pressure with vacancies from the
     * pressure without vacancies for the given lattice density.
     * @param rhoLat lattice density
     * @return vacancy pressure difference
     */
    public double getVacancyPressure(double rhoLat) {
        double c = getVacancyFraction(rhoLat);
        double x = 1 - rhoLat/Math.sqrt(2);
        double x2 = x*x;
        double xp = 1;
        double y2 = 0;
        for (int i=0; i<cdz.length; i++) {
          y2 += cdz[i]*xp;
          xp *= x2;
        }
        return (Math.sqrt(y2) - 3/x)*rhoLat*c;
    }

    /**
     * Returns the Helmholtz free energy difference due to the presence of
     * vacancies at the given lattice density.
     * @param rhoLat
     * @return Helmholtz free energy difference
     */
    public double getVacancyFreeEnergy(double rhoLat) {
        return getVacancyFraction(rhoLat);
    }

    /**
     * Returns the ideal gas free energy for the given density.
     * @param rho Density
     * @return ideal gas free energy
     */
    public static double idFreeEnergy(double rho){return (Math.log(rho)-1.0);}

    public static void main(String[] args) {
        HardSphereSolid hs = new HardSphereSolid();

        if (args.length>0) {
            EOSParams params = new EOSParams();
            ParseArgs.doParseArgs(params, args);
            if (params.density > 0) {
                double rho = params.density;
                double Z = hs.zSolid(rho);
                double A = hs.Asolid(rho, 100000);
                double mu = A + Z;
                double dPdrho = calcdPdrho(rho, params.drho);
                System.out.print(String.format("%20.15f %20.15f %21.14f %21.15f %21.15f\n", rho, Z*rho, A, mu, dPdrho));
            }
            else if (params.pressure > 0) {
                double rho = hs.densityForPressure(params.pressure);
                double Z = hs.zSolid(rho);
                double A = hs.Asolid(rho, 100000);
                double mu = A + Z;
                double dPdrho = calcdPdrho(rho, params.drho);
                System.out.print(String.format("%20.15f %20.15f %21.14f %21.15f %21.15f\n", rho, Z*rho, A, mu, dPdrho));
            }
            else if (params.mu > 0) {
                double rho = hs.densityForMu(params.pressure);
                double Z = hs.zSolid(rho);
                double A = hs.Asolid(rho, 100000);
                double mu = A + Z;
                double dPdrho = calcdPdrho(rho, params.drho);
                System.out.print(String.format("%20.15f %20.15f %21.14f %21.15f %21.15f\n", rho, Z*rho, A, mu, dPdrho));
            }
            return;
        }

        double dr = 0.000001;
        for (int i=100; i<=141; i++) {
            double rho = i*0.01;
            
            double Z = hs.zSolid(rho);
            double A = hs.Asolid(rho, 100000);
            double mu = A + Z;
            System.out.print(String.format("%4.2f %20.15f %21.14f %21.15f %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e\n", rho, Z, A, mu, hs.calcdMudP(rho, dr), (1/rho), (1/rho)*(1 - Z/(hs.calcdPdrho(rho, dr))), hs.calcdmudrho(rho, dr), (1/rho)*(hs.calcdPdrho(rho, dr)), hs.calcdmudrho(rho,dr)/hs.calcdPdrho(rho,dr)));
            
//            System.out.print(String.format("%1.2f %1.2f %20.15f\n", rho1, rho2, -(hs.deltaA(rho1, rho2, 10000)/(1/rho2-1/rho1))));
        }
    }
    
    public static class EOSParams extends ParameterBase {
        public double density = 0;
        public double pressure = 0;
        public double mu = 0;
        public double drho = 0.000001;
    }
}
