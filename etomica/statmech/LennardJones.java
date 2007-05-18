package etomica.statmech;

import etomica.util.Arrays;

/**
 * Implementation of some semi-empirical equations of state for the
 * Lennard-Jones model.
 * 
 * Equations for the fcc solid phase and liquid-fcc saturation conditions are
 * based on the empirical formulas proposed by van der Hoef:
 * 
 * [1] M.A. van der Hoef, "Free energy of the Lennard-Jones solid", Journal of
 * Chemical Physics, vol 113, page 8142 (2000).
 * 
 * Equations for vapor-fcc saturation conditions are based on equations given in:
 * 
 * [2] M.S. van der Hoef, "Gas-solid coexistence of the Lennard-Jones system", 
 * Journal of Chemical Physics, vol 117, page 5093 (2002).
 * 
 * No fluid-phase equations are (yet) implemented here.
 * 
 * This is a library class, containing only static methods, and thus is not
 * intended to be instantiated.
 * 
 * @author David Kofke
 * 
 */
public final class LennardJones {

    private LennardJones() {
        //private constructor prevents instantiation
    }
    /**
     * Molar Helmholtz free energy for the fcc solid, in excess of an ideal gas at the
     * same temperature and density. Total molar free energy is obtained by adding 
	 * T*(ln[rho] - 1) + 3*T*ln[Lambda], where Lambda is the thermal deBroglie
     * wavelength.
     * 
     * More precisely, one should add the term given, plus ln[N!]/N - (ln[N]-1),
     * which goes to zero for N->infinity, where N is the number of atoms
     * (typically specified only for small systems). Note that this detail
     * does not account for all finite-size effects.
     * 
     * Based on Eq. 21 of van der Hoef [1].
     * 
     * @param T
     *            the temperature, in units of epsilon/k (thus input is
     *            kT/epsilon)
     * @param rho
     *            the number density, in units of 1/sigma^3 (thus input is
     *            rho*sigma^3)
     * @return the excess Helmholtz free energy per atom, in units of epsilon 
     */
    public static double aResidualFcc(double T, double rho) {
        double sum = 0.0;
        sum += uStaticFcc(rho) - 1.5 * T * Math.log(T);// harmonic contribution, static-lattice energy plus 3/2 ln(beta)
        sum += T * Uah(T, rho); //anharmonic contribution
        double rhon1 = rho;// rho^(n+1)
        for (int n = 0; n <= 3; n++) {
            sum += T * b[n] / (n + 1) * rhon1;
            rhon1 *= rho;
        }
        sum -= 23.3450759 * T;
        return sum;
    }

    /**
     * Total ensemble-averaged internal potential energy for the fcc solid.
     * Based on Eq. 19 of van der Hoef [1].
     * 
     * @param T
     *            temperature, in units of epsilon/k
     * @param rho
     *            number density, in units of 1/sigma^3
     * @return internal potential energy, in units of epsilon (thus u/epsilon is returned)
     * 
     */
    public static double uFcc(double T, double rho) {
        return uStaticFcc(rho) + 1.5/T + uah(T, rho);
    }
    
    /**
     * Compressibility factor, P/rho/kT for the fcc solid.
     * Based on Eq. 20 of van der Hoef [1].
     * 
     * @param T
     *            temperature, in units of epsilon/k
     * @param rho
     *            number density, in units of 1/sigma^3
     * @return compressibility factor (dimensionless)
     */
    public static double ZFcc(double T, double rho) {
        double sum = 1.0/rho;//ideal-gas value of betaP/rho^2
        sum += rho * (2.0*cStat[0] + 4.0*cStat[1] * rho*rho) / T;// beta dUstatic/drho
        sum += dUahDrho(T, rho);
        double rhon = 1.0;
        for(int n=0; n<=3; n++) {
            sum += b[n] * rhon;
            rhon *= rho;
        }
        return rho*sum;
    }

    /**
     * Potential energy per atom of a static fcc lattice with a LJ particle at
     * each site.
     * 
     * @param rho
     *            number density, in units of 1/sigma^3
     * @return the static-lattice energy, in units of epsilon
     */
    public static double uStaticFcc(double rho) {
        double rho2 = rho * rho;
        return rho2 * (cStat[0] + cStat[1] * rho2);
    }
    
    /**
     * Returns the densities of coexisting (saturated) liquid and fcc phases
     * (i.e., melting and freezing densities). Uses Eqs. (25) and (26) of van
     * der Hoef [1].
     * 
     * @param T
     *            temperature, in units of epsilon/k
     * @return coexistence densities, in units of 1/sigma^3; liquid is element
     *         0, fcc solid is element 1
     */
    public static double[] liquidFccCoexDensities(double T) {
        double beta = 1.0/T;
        double betaM14 = Math.pow(beta, -0.25);
        double rhoLiq = 0.91070 + beta * (-0.25124 + 
                                  beta * (+0.85861 + 
                                  beta * (-1.08918 + 
                                  beta * (+0.63932 + 
                                  beta * (-0.14433)))));
        double rhofcc = 0.92302 + beta * (-0.09218 + 
                                  beta * (+0.62381 + 
                                  beta * (-0.82672 + 
                                  beta * (+0.49124 + 
                                  beta * (-0.10847)))));
        rhoLiq *= betaM14;
        rhofcc *= betaM14;
        return new double[] {rhoLiq, rhofcc};  
    }

    /**
     * Returns the saturation pressure for fluid-fcc coexistence (i.e., the
     * melting or freezing pressure). Uses Eq. (24) of van der Hoef [1], which
     * is based on a formula proposed by Agrawal and Kofke, referenced therein.
     * 
     * @param T
     *            temperature, in units of epsilon/k
     * @return saturation pressure, in units of epsilon/sigma^3 (thus p
     *         sigma^3/epsilon is returned)
     */
    public static double liquidFccCoexPressure(double T) {
        double beta = 1.0/T;
        double betaM54 = Math.pow(beta,-1.25);
        double A = -7.2866;//value recommended by van der Hoef
        double B = -2.9895;//value recommended by van der Hoef
        return betaM54 * Math.exp(-0.4759 * Math.sqrt(beta)) * (16.89 + A * beta + B * beta*beta);
    }
    

    /**
     * Returns the densities of coexisting (saturated) vapor and fcc phases
     * (sublimation densities). Uses Eqs. (2) and (3) of van der Hoef [2].
     * 
     * @param T
     *            temperature, in units of epsilon/k
     * @return coexistence densities, in units of 1/sigma^3; vapor is element
     *         0, fcc solid is element 1
     */
    public static double[] vaporFccCoexDensities(double T) {
        double sum = -cStat[0]*cStat[0]/4.0/cStat[1]/T;
        double Tn = 1.0;
        for(int n=0; n<=5; n++) {
            sum += bSub[n]*Tn;
            Tn *= T;
        }
        double rhoVap = Math.exp(sum);
        sum = Math.sqrt(-cStat[0]/2.0/cStat[1]);
        Tn = T;
        for(int n=1; n<=4; n++) {
            sum += cSub[n-1]*Tn;
            Tn *= T;
        }
        return new double[] {rhoVap, sum};  
    }

    /**
     * Returns the saturation pressure for vapor-fcc coexistence. Uses Eq. (1)
     * of van der Hoef [2].
     * 
     * @param T
     *            temperature, in units of epsilon/k
     * @return saturation pressure, in units of epsilon/sigma^3 (thus p sigma^3/epsilon is returned)
     */
    public static double vaporFccCoexPressure(double T) {
        double sum = -cStat[0]*cStat[0]/4.0/cStat[1]/T;
        double Tn= 1.0;
        for(int n=0; n<=5; n++) {
            sum += aSub[n]*Tn;
            Tn *= T;
        }
        return T * Math.exp(sum);
    }

    /*
     * Implementation of van der Hoef's [1] Eq. 12 for u^{ah}
     */
    private static double uah(double T, double rho) {
        double sum = 0.0;
        double rhon = 1.0;
        for (int n = 0; n <= 2; n++) {
            double betamm = T*T; // beta^(-m) = T^(m), starting at m = 2 is T^2
            for (int m = 2; m <= 5; m++) {
                sum += a[n][m - 2] * rhon * betamm;
                betamm *= T;
            }
            rhon *= rho;
        }
        return -sum;
    }
    
    /*
     * Derivative of Uah wrt rho
     */
    private static double dUahDrho(double T, double rho) {
        double sum = 0.0;
        double rhon1 = 1.0;
        for (int n = 1; n <= 2; n++) {
            double beta1m = T; // beta^(1-m) = T^(m-1), starting at m = 2 is T
            for (int m = 2; m <= 5; m++) {
                sum += n * a[n][m - 2] / (m - 1) * rhon1 * beta1m;
                beta1m *= T;
            }
            rhon1 *= rho;
        }
        return -sum;
    }

    /*
     * Implementation of van der Hoef's [1] Eq. (14) for U^{ah}
     */
    private static double Uah(double T, double rho) {
        double sum = 0.0;
        double rhon = 1.0;
        for (int n = 0; n <= 2; n++) {
            double beta1m = T; // beta^(1-m) = T^(m-1), starting at m = 2 is T
            for (int m = 2; m <= 5; m++) {
                sum += a[n][m - 2] / (m - 1) * rhon * beta1m;
                beta1m *= T;
            }
            rhon *= rho;
        }
        return -sum;
    }
    
    public static void main(String[] args) {
        double T = 1.0;
        double rho = 0.962;
        double Tm = 0.689;

        double pSat = liquidFccCoexPressure(T);
        double[] rhoSat = liquidFccCoexDensities(T);
        System.out.println("Input T, rho: "+ T + "  "+rho);
        System.out.println("Residual Helmholtz: " + aResidualFcc(T, rho));
        double betaA = (aResidualFcc(T, rho) + (-1 + Math.log(rho)));
        System.out.println("betaA: " + betaA);
        System.out.println("betaA anharmonic: " + Uah(T, rho));
        System.out.println("Potential energy: " + uFcc(T, rho));
        System.out.println("Lattice energy: " + uStaticFcc(rho));
        System.out.println("Compressibility factor: " + ZFcc(T, rho));
        System.out.println("Liquid-fcc coexistence pressure: " + pSat);
        System.out.println("liquid, fcc coexistence densities: "+ Arrays.toString(rhoSat));
        double rhs = aResidualFcc(T,rho)/T - uStaticFcc(rho)/T + 1.5*Math.log(T) - Uah(T,rho) 
                    + 1.5*Math.log(2*Math.PI);
        System.out.println("rhs of Eq. 23: " + rhs);
        double zfcc = pSat/rhoSat[1]/T;
        System.out.println("Two routes to Z: "+zfcc+" "+ZFcc(T,rhoSat[1]));
        double psub = vaporFccCoexPressure(T);
        System.out.println("Sublimation pressure: " + psub);
        double[] rhosub = vaporFccCoexDensities(T);
        System.out.println("vapor, fcc coexistence densities: " + Arrays.toString(rhosub));
        
//        for(double Tr=0.2; Tr<=1.01; Tr+=0.2) {
            for(double Tr=0.0; Tr<=1.01; Tr+=0.05) {
//            for(rho=0.9; rho<=1.7; rho+=0.1) {
                T = Tr*Tm;
                double a = (aResidualFcc(T, rho) + T*(-1 + Math.log(rho)));
//                System.out.println(Tr + "\t" + a);
                System.out.println(a);
//                System.out.println("T, rho, Uah: "+ T + "\t" + rho + "\t" + Uah(T,rho));
//            }
        }
    }
    
    //Constants for van der Hoef's formulas, ref[1]
    private static final double[][] a = new double[][] {
        { -8.2151768, 12.070686, -6.6594615, 1.3211582 },
        { 13.404069, -20.632066, 11.564825, -2.3064801 },
        { -5.5481261, 8.8465978, -5.0258631, 1.0070066 } };
    private static final double[] b = new double[] { 69.833875, -132.86963,
        97.438593, -25.848057 };
    private static final double[] cStat = new double[] {-14.45392093, 6.065940096};

    //Constants for van der Hoef's formulas, ref[2]
    private static final double[] aSub = new double[] {10.48631, -16.81771, 35.10031, -48.76487, 36.39136, -10.90788};
    private static final double[] bSub = new double[] {10.49107, -16.83818, 34.95754, -47.66412, 33.90821, -9.011857};
    private static final double[] cSub = new double[] {-0.134343, -0.0950795, 0.137215, -0.161890};


}
