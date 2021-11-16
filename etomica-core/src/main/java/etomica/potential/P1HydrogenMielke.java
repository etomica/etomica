/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomHydrogen;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.BohrRadius;
import etomica.units.Hartree;

public class P1HydrogenMielke {
//    1 Eh (hartree) = 27.2113961eV = 627.5096 kcal/mol = 219474.7 cm^-1 

    protected final static double c6 = -6.499027;
    protected final static double c8 = -124.3991;
    protected final static double c10 = -3285.828;
    protected final static double beta = 0.2;
    protected final static double re = 1.401;
    protected final static double alpha = 3.980917850296971;
    protected final static double[] xlp = {-1.709663318869429E-1,-6.675286769482015E-1,-1.101876055072129E+0,-1.106460095658282E+0,-7.414724284525231E-1,-3.487882731923523E-1,-1.276736255756598E-1,-5.875965709151867E-2,-4.030128017933840E-2,-2.038653237221007E-2,-7.639198558462706E-4,2.912954920483885E-3,-2.628273116815280E-4,-4.622088855684211E-4,1.278468948126147E-4,-1.157434070240206E-5,-2.609691840882097E-12};

    public P1HydrogenMielke(Space space) {
    }

    public double u(double rSim) {       
        // H2 singlet potential for FCI/CBS
    	if (rSim < 0) return Double.POSITIVE_INFINITY;
        double rBohr = BohrRadius.UNIT.fromSim(rSim);
        double damp=1.0-Math.exp(-beta*rBohr*rBohr);
        double inverseRm = damp/rBohr ;
        double inverseRm2 = inverseRm*inverseRm;
        double inverseRm6 = inverseRm2*inverseRm2*inverseRm2;
        double inverseRm8 = inverseRm6*inverseRm2;
        double inverseRm10 = inverseRm8*inverseRm2;
        double prefac = Math.exp(- alpha*(rBohr-re));
        double uLong = c6*inverseRm6 + c8*inverseRm8 + c10*inverseRm10;
        double uShort = 0.00;        
        for (int i=0; i<17; i++) {
            uShort += xlp[i]*prefac;
            prefac *= (rBohr-re);
        }
        double uTotal = Hartree.UNIT.toSim(uLong + uShort);
        return uTotal;
    }
    
    public double du (double rSim) {// returns rSim*dudr
    	double rBohr = BohrRadius.UNIT.fromSim(rSim);
    	double r2 = rBohr*rBohr;
    	double z = beta*r2;
        double y = Math.exp(-z);
        double inverseRm = (1.0 - y)/rBohr ;
        double inverseRm2 = inverseRm*inverseRm;
        double inverseRm5 = inverseRm2*inverseRm2*inverseRm;
        double inverseRm7 = inverseRm5*inverseRm2;
        double inverseRm9 = inverseRm7*inverseRm2;        
        double dinverseRmdr = (2*z*y + y - 1)/r2;
        double duLongdr = dinverseRmdr*(6*c6*inverseRm5 + 8*c8*inverseRm7 + 10*c10*inverseRm9);
        double duShortdr = -alpha*xlp[0];
        double fac = 1.0;
        for (int i=1; i<17; i++) {
            duShortdr += xlp[i]*fac*(i - alpha*(rBohr-re));
            fac *= (rBohr-re);
        }
        duShortdr *= Math.exp(-alpha*(rBohr-re));
        double duTotal = Hartree.UNIT.toSim(duLongdr + duShortdr)/BohrRadius.UNIT.toSim(1);         
        return rSim*duTotal;
    }   
    
    public double d2u(double rSim) {// returns rSim^2*d2udr2
        double rBohr = BohrRadius.UNIT.fromSim(rSim);
        double r2 = rBohr*rBohr;
    	double z = beta*r2;
        double y = Math.exp(-z);        
        double inverseRm = (1.0 - y)/rBohr ;
        double inverseRm2 = inverseRm*inverseRm;
        double inverseRm4 = inverseRm2*inverseRm2;
        double inverseRm5 = inverseRm4*inverseRm;
        double inverseRm6 = inverseRm5*inverseRm;
        double inverseRm7 = inverseRm6*inverseRm;
        double inverseRm8 = inverseRm7*inverseRm;
        double inverseRm9 = inverseRm8*inverseRm;
        double dinverseRmdr = (2*z*y + y - 1)/r2;
        double d2inverseRmdr2 = (2 - 2*y - 2*z*y - 4*z*z*y)/(r2*rBohr);        
        double d2uL = d2inverseRmdr2*(6*c6*inverseRm5 + 8*c8*inverseRm7 + 10*c10*inverseRm9) + dinverseRmdr*dinverseRmdr*(30*c6*inverseRm4 + 56*c8*inverseRm6 + 90*c10*inverseRm8);
        double d2uS = alpha*alpha*xlp[0] + (alpha*alpha*(rBohr-re) - 2*alpha)*xlp[1];
        double fac = 1;
        for (int i=2; i<17; i++) {
            d2uS += xlp[i]*fac*(alpha*alpha*(rBohr-re)*(rBohr-re) + i*(i-1) - 2*i*alpha*(rBohr-re));
            fac *= (rBohr-re);
        }
        d2uS *= Math.exp(-alpha*(rBohr-re));
        double d2uT = Hartree.UNIT.toSim(d2uL + d2uS)/(BohrRadius.UNIT.toSim(1)*BohrRadius.UNIT.toSim(1));
        return rSim*rSim*d2uT;
    }
    
    
    public static class P1HydrogenMielkeAtomic extends P1HydrogenMielke implements IPotentialField {
        public P1HydrogenMielkeAtomic(Space space) {
            super(space);     
        }

        public double getRange() {
            return Double.POSITIVE_INFINITY;
        }

        @Override
        public double u(IAtom atom) {
            double bL = ((AtomHydrogen)atom).getBondLength();
            return u(bL);
        }

        @Override
        public double udu(IAtom atom, Vector f) {
            return u(atom);
        }
    }
    
    public static class P2HydrogenMielkeAtomic extends P1HydrogenMielke implements IPotentialAtomic {
    	public P2HydrogenMielkeAtomic(Space space) {
    		super(space);
    	}

        public double getRange() {
            return Double.POSITIVE_INFINITY;
        }

    	public double energy(IAtomList atoms) {
    		IAtom a0 = atoms.get(0);
    		IAtom a1 = atoms.get(1);
            double bL = Math.sqrt(a0.getPosition().Mv1Squared(a1.getPosition()));
            double f = u(bL);
            return f;
    	}
    }
    
}
