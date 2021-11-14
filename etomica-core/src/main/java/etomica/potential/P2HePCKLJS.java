/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.util.Constants;

/**
 *
 * Ab initio pair potential for helium from Przybytek et al. (2010) Phys. Rev. Lett. 104, 183003.   
 *
 * Potential is speciously negative at separations less than 0.3 a0.
 * Second derivative (used for quadratic Feymann-Hibbs potential) goes through maximum at 0.2A~=0.4a0. 
 * We apply a hard core of 0.4 a0.
 *
 * @author Kate Shaul and Andrew Schultz
 */
public class P2HePCKLJS extends Potential2SoftSpherical {

    public static Potential2Soft makeTruncated(Space space, double sigma, TruncationFactory tf) {
        return tf.make(new P2HePCKLJS(space, sigma));
    }

    public P2HePCKLJS(Space space) {
        this(space, 0);
    }

    public P2HePCKLJS(Space space, double sigma) {
        super(space);

        double W4 = 0.35322e-04 / (alpha * alpha);
        double AS3 = 0.577235e-06 / (alpha * alpha * alpha);
        double polarizability = 1.38319217440;

        double K7 = 23.0/(4.0*Math.PI)*polarizability*polarizability/alpha;
        double ratio = alpha*K7/C6BO;

        B = new double[7];
        B[1] = 8.454943177941253;
    	B[2] = 16.006586066260556;
    	B[3] = 10.378373954734820;
    	B[4] = 3.515803817223855;
    	B[5] = 0.591502377533792;
    	B[6] = 0.059455768329599;

        A = new double[7];
        A[1] = B[1];
        A[2] = B[2]-W4/C6BO;
        A[3] = B[3]-B[1]*W4/C6BO+AS3/C6BO;
        A[4] = ratio*B[5];
        A[5] = ratio*B[6];
        A[6] = 0;

        C = new double[17];
        C[3] = 0.000000577235; 
    	C[4] = -0.000035322; 
    	C[5] = 0.000001377841; 
    	C[6] = 1.461830; 
    	C[8] = 14.12350; 
    	C[10] = 183.7497; 
    	C[11] = -0.7674*100; 
    	C[12] = 0.3372*1e4; 
    	C[13] = -0.3806*1e4; 
    	C[14] = 0.8534*1e5; 
    	C[15] = -0.1707*1e6;  
    	C[16] = 0.286*1e7; 
   
    	P = new double[] { -25.4701669416621, 269.244425630616, -56.3879970402079};
    	Q = new double[] {38.7957487310071, -2.76577136772754};

    	aErr = new double[]{2.456, 1.100, 0.4381};
    	cErr = new double[]{0.16702e-3, 0.4524e-5, 0.1843e-7};

    	setErrMulti(sigma);
    }

    public void setErrMulti(double newErrMult) {
    	errMult = newErrMult;
    }

    /**
     * The energy u.
     */
    public double u(double r2) {

    	double r = Math.sqrt(r2);
    	r = BohrRadius.UNIT.fromSim(r);

    	if (r < sigmaHC) {
    		return Double.POSITIVE_INFINITY;
    	}

    	double u1 = (P[0] + P[1]*r + P[2]*r*r)*Math.exp(-a*r);

    	double u2 = (Q[0] + Q[1]*r)*Math.exp(-b*r);
    	
    	double invr = 1.0/r;
    	
    	double br = eta*r;
        double m = Math.exp(-br);
        double term = 1.0;
        double sum = term;
        double invri = invr;
        double u3 = 0;
        for (int i=1; i<17; i++) {
            term *= br/i;
            sum = sum + term;
            u3 += (-1.0+m*sum)*C[i]*invri;
            invri *= invr;
        }

    	/// damp_ret ////

        double sumA = 1.0;
        double sumB = 1.0;

        double x = alpha*r;
        double xn=1.0;
        for (int n=1;n<7;n++) {
            xn *= x;
            sumA += A[n]*xn;
            sumB += B[n]*xn;
        }
        double g = sumA/sumB;
        
        double invr3 = invr*invr*invr;
    	double Vret = (C[3] + C[4]*invr + C6BO*(1.0-g)*invr3)*invr3;
    	double uSigma = 0;
    	if (errMult != 0) {
    		double s = 0;
        	for (int i=0; i<3; i++) {
        		s += cErr[i]*Math.exp(-aErr[i]*r);
        	}
    		uSigma = errMult*s;
    	}
 
    	return  Hartree.UNIT.toSim(u1 + u2 + u3 + Vret + uSigma);
    	
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) {
        
    	double r = Math.sqrt(r2);
    	r = BohrRadius.UNIT.fromSim(r);

    	if (r < sigmaHC) {
    		return 0;
    	}

    	double u1 = (P[0] + P[1]*r + P[2]*r*r)*Math.exp(-a*r);
    	double du1dr = (P[1] + P[2]*2.0*r)*Math.exp(-a*r) - a*u1;
    	double u2 = (Q[0] + Q[1]*r)*Math.exp(-b*r);
    	double du2dr = Q[1]*Math.exp(-b*r) - b*u2;

    	double invr = 1.0/r;

    	double br = eta*r;
        double m = Math.exp(-br);
        double dmdr = -eta*m;
        double term = 1.0;
        double sum = term;
        double dsumdr = 0;
    	double du3dr = 0;
        double invri = invr;
        for (int i=1; i<17; i++) {
            term *= br/i;
            sum = sum + term;
            dsumdr = dsumdr + term*invr*i;
            du3dr += (dmdr*sum + m*dsumdr)*C[i]*invri + (-1.0+m*sum)*C[i]*invri*invr*(-i);
            invri *= invr;
        }


    	/// damp_ret ////

        double sumA = 1.0;
        double sumB = 1.0;
        double dsumAdr = 0;
        double dsumBdr = 0;
        double x = alpha*r;
        double xn=1.0;
        for (int n=1;n<7;n++) {
            xn *= x;
            sumA += A[n]*xn;
            sumB += B[n]*xn;
            dsumAdr += A[n]*xn*invr*n;
            dsumBdr += B[n]*xn*invr*n;
        }
        double g = sumA/sumB;
        double dgdr = dsumAdr/sumB - sumA/sumB/sumB*dsumBdr;
        double invr3 = invr*invr*invr;
    	//double Vret = (C[3] + C[4]*invr + C6BO*(1.0-g)*invr3)*invr3;
    	double dVretdr = (-3.0*C[3] -4.0*C[4]*invr -6.0*C6BO*(1.0-g)*invr3)*invr3*invr + C6BO*(-dgdr)*invr3*invr3;
        
    	double duSigma = 0;
    	if (errMult != 0) {
    		double s = 0;
        	for (int i=0; i<3; i++) {
        		s += -aErr[i]*cErr[i]*Math.exp(-aErr[i]*r);
        	}
    		duSigma = errMult*s;
    	}

    	double dudr = Hartree.UNIT.toSim(du1dr + du2dr + du3dr + dVretdr + duSigma);
    	
        return r*dudr;
        
        
    
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        double r = Math.sqrt(r2);
        r = BohrRadius.UNIT.fromSim(r);

        if (r < sigmaHC) {
            return 0;
        }

        double expar = Math.exp(-a*r);
        double u1 = (P[0] + P[1]*r + P[2]*r*r)*expar;
//        double du1dr = (P[1] + P[2]*2.0*r)*Math.exp(-a*r) - a*u1;
        double d2u1dr2 = 2*P[2]*expar -2*a*(P[1] + P[2]*2.0*r)*expar + a*a*u1;
        double expbr = Math.exp(-b*r);
        double u2 = (Q[0] + Q[1]*r)*expbr;
//        double du2dr = Q[1]*Math.exp(-b*r) - b*(Q[0] + Q[1]*r)*Math.exp(-b*r);
        double d2u2dr2 = -2*b*Q[1]*expbr + b*b*u2;
        
        double invr = 1.0/r;
        invr = 1.0/r;

        double m = Math.exp(-eta*r);
        double dmdr = -eta*m;
        double d2mdr2 = -eta*dmdr;
        double term = 1.0;
        double sum = term;
        double dsumdr = 0;
        double d2sumdr2 = 0;
        double d2u3dr2 = 0;
        double invri = invr;
        for (int i=1; i<17; i++) {
            term *= eta*r/i;
            sum = sum + term;
            dsumdr = dsumdr + term*invr*i;
            d2sumdr2 += term*invr*invr*i*(i-1);
//            u3 += (-1.0+m*sum)*C[i]*invri;
//            du3dr += (dmdr*sum + m*dsumdr)*C[i]*invri + (-1.0+m*sum)*C[i]*invri*invr*(-i);
            d2u3dr2 += (d2mdr2*sum + 2*dmdr*dsumdr + m*d2sumdr2)*C[i]*invri
                     + 2*(dmdr*sum + m*dsumdr)*C[i]*invri*invr*(-i)
                     + (-1.0+m*sum)*C[i]*invri*invr*invr*(-i)*(-i-1);
            invri *= invr;
        }

        /// damp_ret ////

        double sumA = 1.0;
        double sumB = 1.0;
        double dsumAdr = 0;
        double dsumBdr = 0;
        double d2sumAdr2 = 0;
        double d2sumBdr2 = 0;
        double x = alpha*r;
        double xn=1.0;
        for (int n=1;n<7;n++) {
            xn *= x;
            double termA = A[n]*xn;
            double termB = B[n]*xn;
            sumA += termA;
            sumB += termB;
            termA *= invr*n;
            termB *= invr*n;
            dsumAdr += termA;
            dsumBdr += termB;
            termA *= invr*(n-1);
            termB *= invr*(n-1);
            d2sumAdr2 += termA;
            d2sumBdr2 += termB;
        }
        double g = sumA/sumB;
        double dgdr = dsumAdr/sumB - sumA/sumB/sumB*dsumBdr;
        double d2gdr2 = (d2sumAdr2 - (2*dsumAdr*dsumBdr + sumA*d2sumBdr2 - 2*sumA*dsumBdr*dsumBdr/sumB)/sumB)/sumB;
        double invr3 = invr*invr*invr;
//        double Vret = (C[3] + C[4]*invr + C6BO*(1.0-g)*invr3)*invr3;
//        double dVretdr = (-3.0*C[3] -4.0*C[4]*invr -6.0*C6BO*(1.0-g)*invr3)*invr3*invr + C6BO*(-dgdr)*invr3*invr3;
        double d2Vretdr2 = (12.0*C[3] +20.0*C[4]*invr +42.0*C6BO*(1.0-g)*invr3)*invr3*invr*invr + C6BO*(-d2gdr2 + 12*dgdr*invr)*invr3*invr3;

    	double d2uSigma = 0;
    	if (errMult != 0) {
    		double s = 0;
        	for (int i=0; i<3; i++) {
        		s += aErr[i]*aErr[i]*cErr[i]*Math.exp(-aErr[i]*r);
        	}
    		d2uSigma = errMult*s;
    	}

        double d2udr2 = Hartree.UNIT.toSim(d2u1dr2 + d2u2dr2 + d2u3dr2 + d2Vretdr2 + d2uSigma);
        return r*r*d2udr2;
    }
            
    /**
     * Integral used for corrections to potential truncation.
     */
    public double integral(double rC) {

        return 0;  //complete LRC is obtained by multiplying by N1*N2/V
    }
    
    public P2HeQFH makeQFH(double temperature) {
        return new P2HeQFH(temperature);
    }    
    
    /**
     * This inner class can calculates the Feynman-Hibbs semiclassical
     * approximation for the potential.  Results should be the same as using
     * P2EffectiveFeynmanHibbs, but should be 2x faster, because u, du and d2u
     * are not called.  Much of the work is duplicated between those methods,
     * and they are all computed at the same time in this class.
     *
     * @author Andrew Schultz
     */
    public class P2HeQFH implements Potential2Spherical {

        protected final double temperature;
        protected final double mass = 4.002602;
        protected double fac;

        public P2HeQFH(double temperature) {
            this.temperature = temperature;
            double hbar = Constants.PLANCK_H/(2*Math.PI);
            fac = hbar*hbar/(24*mass/2)/temperature;
        }

        public double energy(IAtomList atoms) {
            dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
            boundary.nearestImage(dr);
            return u(dr.squared());
        }

        public double getRange() {
            return P2HePCKLJS.this.getRange();
        }

        public double u(double r2) {
            double r = Math.sqrt(r2);
            r = BohrRadius.UNIT.fromSim(r);

            if (r < sigmaHC) {
                return Double.POSITIVE_INFINITY;
            }

            double expar = Math.exp(-a*r);
            double expbr = Math.exp(-b*r);
            
            double br = eta*r;
            double m = Math.exp(-br);

            double u1 = (P[0] + P[1]*r + P[2]*r*r)*expar;
            double du1dr = (P[1] + P[2]*2.0*r)*expar - a*u1;
            double d2u1dr2 = 2*P[2]*expar -2*a*(P[1] + P[2]*2.0*r)*expar + a*a*u1;
            double u2 = (Q[0] + Q[1]*r)*expbr;
            double du2dr = Q[1]*expbr - b*(Q[0] + Q[1]*r)*expbr;
            double d2u2dr2 = -2*b*Q[1]*expbr + b*b*u2;
          
            double invr = 1.0/r;

            double u3 = 0;
            double du3dr = 0;
            double dmdr = -eta*m;
            double d2mdr2 = -eta*dmdr;
            double term = 1.0;
            double sum = term;
            double dsumdr = 0;
            double d2sumdr2 = 0;
            double d2u3dr2 = 0;
            double invri = invr;
            for (int i=1; i<17; i++) {
                term *= br/i;
                sum = sum + term;
                dsumdr = dsumdr + term*invr*i;
                d2sumdr2 += term*invr*invr*i*(i-1);
                u3 += (-1.0+m*sum)*C[i]*invri;
                du3dr += (dmdr*sum + m*dsumdr)*C[i]*invri + (-1.0+m*sum)*C[i]*invri*invr*(-i);
                d2u3dr2 += (d2mdr2*sum + 2*dmdr*dsumdr + m*d2sumdr2)*C[i]*invri
                       + 2*(dmdr*sum + m*dsumdr)*C[i]*invri*invr*(-i)
                       + (-1.0+m*sum)*C[i]*invri*invr*invr*(-i)*(-i-1);
                invri *= invr;
            }

            /// damp_ret ////

            double sumA = 1.0;
            double sumB = 1.0;
            double dsumAdr = 0;
            double dsumBdr = 0;
            double d2sumAdr2 = 0;
            double d2sumBdr2 = 0;
            double x = alpha*r;
            double xn=1.0;
            for (int n=1;n<7;n++) {
                xn *= x;
                double termA = A[n]*xn;
                double termB = B[n]*xn;
                sumA += termA;
                sumB += termB;
                termA *= invr*n;
                termB *= invr*n;
                dsumAdr += termA;
                dsumBdr += termB;
                termA *= invr*(n-1);
                termB *= invr*(n-1);
                d2sumAdr2 += termA;
                d2sumBdr2 += termB;
            }
            double g = sumA/sumB;
            double dgdr = dsumAdr/sumB - sumA/sumB/sumB*dsumBdr;
            double d2gdr2 = (d2sumAdr2 - (2*dsumAdr*dsumBdr + sumA*d2sumBdr2 - 2*sumA*dsumBdr*dsumBdr/sumB)/sumB)/sumB;
            double invr3 = invr*invr*invr;
            double Vret = (C[3] + C[4]*invr + C6BO*(1.0-g)*invr3)*invr3;
            double dVretdr = (-3.0*C[3] -4.0*C[4]*invr -6.0*C6BO*(1.0-g)*invr3)*invr3*invr + C6BO*(-dgdr)*invr3*invr3;
            double d2Vretdr2 = (12.0*C[3] +20.0*C[4]*invr +42.0*C6BO*(1.0-g)*invr3)*invr3*invr*invr + C6BO*(-d2gdr2 + 12*dgdr*invr)*invr3*invr3;

        	double uSigma = 0, duSigma = 0, d2uSigma = 0;
        	if (errMult != 0) {
        		double s = 0, ds = 0, d2s = 0;
            	for (int i=0; i<3; i++) {
            		double t = cErr[i]*Math.exp(-aErr[i]*r);
            		s += t;
            		t *= -aErr[i];
            		ds += t;
            		t *= -aErr[i];
            		d2s += t;
            	}
        		uSigma = errMult*s;
        		duSigma = errMult*ds;
        		d2uSigma = errMult*d2s;
        	}

            double uc = Hartree.UNIT.toSim(u1 + u2 + u3 + Vret + uSigma);
            double duc = Hartree.UNIT.toSim(r*(du1dr + du2dr + du3dr + dVretdr + duSigma));
            double d2uc = Hartree.UNIT.toSim(r*r*(d2u1dr2 + d2u2dr2 + d2u3dr2 + d2Vretdr2 + d2uSigma));

            if (uc == Double.POSITIVE_INFINITY || d2uc == Double.POSITIVE_INFINITY) { return Double.POSITIVE_INFINITY; }
            double u = uc + (fac/r2)*(d2uc + 2*duc);

            // if the classical potential is repulsive, the semiclassical potential
            // should be more repulsive.  In nearly all cases, it is, but this often
            // fails for very short separations.  Just enforce it here.
            if (uc > 0 && u < uc) return uc;
            return u;
        }
    }
    
    public P2HeTI makeTI(double temperature) {
        return new P2HeTI(temperature);
    }
    
    public class P2HeTI implements Potential2Spherical {

        protected final double temperature;
        protected final double mass = 4.002602;
        protected final double fac;

        public P2HeTI(double temperature) {
            this.temperature = temperature;
            double hbar = Constants.PLANCK_H/(2*Math.PI);
            fac = hbar*hbar/(24*mass/2)/(temperature*temperature);
        }

        public double energy(IAtomList atoms) {
            dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
            boundary.nearestImage(dr);
            return u(dr.squared());
        }

        public double getRange() {
            return P2HePCKLJS.this.getRange();
        }

        public double u(double r2) {
            double r = Math.sqrt(r2);            
            r = BohrRadius.UNIT.fromSim(r);

            if (r < sigmaHC) {
                return Double.POSITIVE_INFINITY;
            }

            double expar = Math.exp(-a*r);
            double expbr = Math.exp(-b*r);
            
            double br = eta*r;
            double m = Math.exp(-br);

            double u1 = (P[0] + P[1]*r + P[2]*r*r)*expar;
            double du1dr = (P[1] + P[2]*2.0*r)*expar - a*u1;
            double d2u1dr2 = 2*P[2]*expar -2*a*(P[1] + P[2]*2.0*r)*expar + a*a*u1;
            double u2 = (Q[0] + Q[1]*r)*expbr;
            double du2dr = Q[1]*expbr - b*(Q[0] + Q[1]*r)*expbr;
            double d2u2dr2 = -2*b*Q[1]*expbr + b*b*u2;
          
            double invr = 1.0/r;

            double u3 = 0;
            double du3dr = 0;
            double dmdr = -eta*m;
            double d2mdr2 = -eta*dmdr;
            double term = 1.0;
            double sum = term;
            double dsumdr = 0;
            double d2sumdr2 = 0;
            double d2u3dr2 = 0;
            double invri = invr;
            for (int i=1; i<17; i++) {
                term *= br/i;
                sum = sum + term;
                dsumdr = dsumdr + term*invr*i;
                d2sumdr2 += term*invr*invr*i*(i-1);
                u3 += (-1.0+m*sum)*C[i]*invri;
                du3dr += (dmdr*sum + m*dsumdr)*C[i]*invri + (-1.0+m*sum)*C[i]*invri*invr*(-i);
                d2u3dr2 += (d2mdr2*sum + 2*dmdr*dsumdr + m*d2sumdr2)*C[i]*invri
                       + 2*(dmdr*sum + m*dsumdr)*C[i]*invri*invr*(-i)
                       + (-1.0+m*sum)*C[i]*invri*invr*invr*(-i)*(-i-1);
                invri *= invr;
            }

            /// damp_ret ////

            double sumA = 1.0;
            double sumB = 1.0;
            double dsumAdr = 0;
            double dsumBdr = 0;
            double d2sumAdr2 = 0;
            double d2sumBdr2 = 0;
            double x = alpha*r;
            double xn=1.0;
            for (int n=1;n<7;n++) {
                xn *= x;
                double termA = A[n]*xn;
                double termB = B[n]*xn;
                sumA += termA;
                sumB += termB;
                termA *= invr*n;
                termB *= invr*n;
                dsumAdr += termA;
                dsumBdr += termB;
                termA *= invr*(n-1);
                termB *= invr*(n-1);
                d2sumAdr2 += termA;
                d2sumBdr2 += termB;
            }
            double g = sumA/sumB;
            double dgdr = dsumAdr/sumB - sumA/sumB/sumB*dsumBdr;
            double d2gdr2 = (d2sumAdr2 - (2*dsumAdr*dsumBdr + sumA*d2sumBdr2 - 2*sumA*dsumBdr*dsumBdr/sumB)/sumB)/sumB;
            double invr3 = invr*invr*invr;
            double Vret = (C[3] + C[4]*invr + C6BO*(1.0-g)*invr3)*invr3;
            double dVretdr = (-3.0*C[3] -4.0*C[4]*invr -6.0*C6BO*(1.0-g)*invr3)*invr3*invr + C6BO*(-dgdr)*invr3*invr3;
            double d2Vretdr2 = (12.0*C[3] +20.0*C[4]*invr +42.0*C6BO*(1.0-g)*invr3)*invr3*invr*invr + C6BO*(-d2gdr2 + 12*dgdr*invr)*invr3*invr3;

            double uc = Hartree.UNIT.toSim(u1 + u2 + u3 + Vret);
            double duc = Hartree.UNIT.toSim(r*(du1dr + du2dr + du3dr + dVretdr));
            double d2uc = Hartree.UNIT.toSim(r*r*(d2u1dr2 + d2u2dr2 + d2u3dr2 + d2Vretdr2));

            if (uc == Double.POSITIVE_INFINITY || d2uc == Double.POSITIVE_INFINITY) { return Double.POSITIVE_INFINITY; }
            double u = uc + (fac/r2)*duc*duc;
            
//            System.out.print(ry+" "+u+" ");
            return u;
        }
    }
    
    public static void main(String[] args) {

    	Space space = Space3D.getInstance();
    	P2HePCKLJS p2 = new P2HePCKLJS(space);

    	//Test separations in fortran code, in Hartrees
    	double[] VTotRet = new double [] {9.07505678537836e-01,3.30439541259472e-01,1.14486105544985e-01,3.78887051832271e-02,1.19356572795740e-02,3.51858623016233e-03,9.26991759357220e-04,1.85096029792210e-04,-1.44666519197000e-06,-2.89977763063000e-05,-3.48209516983400e-05,-3.06464226675000e-05,-2.18216351894100e-05,-1.46331763646000e-05,-6.54141208063000e-06,-3.13153787265000e-06,-5.24252133480000e-07};
    	double[] rTest = new double [] {1.00 , 1.50 , 2.00 , 2.50 , 3.00 , 3.50 , 4.00 , 4.50 , 5.00 , 5.30 , 5.60 , 6.00 , 6.50 , 7.00 , 8.00 , 9.00 , 12.00}; 	
    	
    	double r;double u;
     	System.out.println();
    	System.out.println("r(a0) \t V(Tot+Ret) testEh.f90 \t V(Tot+Ret) here \t here-testEh.f90 (K)   ");
     	for (int i=0;i<rTest.length;i++) {
    		r = BohrRadius.UNIT.toSim(rTest[i]); //Angstrom
    		u = Hartree.UNIT.fromSim(p2.u(r*r)); // sim
    		double V = VTotRet[i];
    	 
    		System.out.println(rTest[i]+"   \t"+V+"  \t" +u+"  \t" +(u-V));
    	}
     	
     	r=5;
     	while (r<6) {

     		r = r + 0.1;
    		double rA = BohrRadius.UNIT.toSim(r); //Angstrom
    		u = Hartree.UNIT.fromSim(p2.u(rA*rA)); // sim
    	 
    		System.out.println(r+"   \t"+u);
    	}

/*     	long t1 = System.currentTimeMillis();
        for (r = 1; r<10; r+=0.0000001) {
            p2.u(r*r);
        }
        long t2 = System.currentTimeMillis();
        System.out.println((t2-t1)/1000.0);*/
        
     	//Sanity check on rdu/dr
     	System.out.println();
     	r=5; 
     	double rmin = r;
     	double u1 = p2.u(r*r);
     	double delr = 0.000001;
     	System.out.println("r(A) \t\t rdudr \t\t\t rdeludelr \t\t\t\t u2");
     	while (r<(rmin+0.00001)) {

     		r = r + delr;
    		double u2 = p2.u(r*r); //u2 expects Angstroms - r is Angstroms
    		double rdudr = p2.du(r*r);
    		double rdeludelr = r*(u2-u1)/delr;
    		
    		//System.out.println(r+"   \t"+dudr+"  \t" +u2);
    		System.out.println(r+"   \t"+rdudr+"  \t" +rdeludelr+"  \t" +u2);
    		u1=u2;
    	}
     	
     	//Sanity check on r2d2u/dr2
     	System.out.println();
     	r = 2; rmin=r; 

     	double dudr1 = p2.du(r*r)/r;
     	delr = 0.000001;
     	System.out.println("r(A) \t\t d2udr2 \t\t\t del2udelr2 \t\t\t\t u2");
     	while (r<(rmin+0.00001)) {

     		r = r + delr;
    		double dudr2 = p2.du(r*r)/r;
    		double d2udr2 = p2.d2u(r*r)/r/r;
    		double del2udelr2 = (dudr2-dudr1)/delr;
    		
    		//System.out.println(r+"   \t"+dudr+"  \t" +u2);
    		System.out.println(r+"   \t"+d2udr2+"  \t" +del2udelr2+"  \t");
    		dudr1=dudr2;
    	}
    }
    
    private static final long serialVersionUID = 1L;
    protected final double[] C, A, B, P, Q;
    protected static final double C6BO = 1.460977837725; //pulled from potentials.f90...
    protected static final double alpha = 1.0/137.036; //fsalpha from potentials.f90...
    protected static final double a = 3.64890303652830;
    protected static final double b = 2.36824871743591;
    protected static final double eta = 4.09423805117871;
    protected static double sigmaHC = 0.4; //bohr radii
    protected final double[] aErr, cErr;
    protected double errMult;
}
