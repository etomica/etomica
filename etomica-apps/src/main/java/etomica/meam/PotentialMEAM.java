/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.PotentialN;
import etomica.potential.PotentialSoft;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.util.Arrays;

/**
 * @author Kate Schadel
 */
public class PotentialMEAM extends PotentialN implements PotentialSoft {

    public static final int RHOj0 = 0;
    public static final int RHOj1x = 1;
    public static final int RHOj1y = 2;
    public static final int RHOj1z = 3;
    public static final int RHOj2xx = 4;
    public static final int RHOj2xy = 5;
    public static final int RHOj2xz = 6;
    public static final int RHOj2yy = 7;
    public static final int RHOj2yz = 8;
    public static final int RHOj2zz = 9;
    public static final int RHOj2 = 10;
    public static final int RHOj3xxx = 11;
    public static final int RHOj3xxy = 12;
    public static final int RHOj3xxz = 13;
    public static final int RHOj3xyy = 14;
    public static final int RHOj3xyz = 15;
    public static final int RHOj3xzz = 16;
    public static final int RHOj3yyy = 17;
    public static final int RHOj3yyz = 18;
    public static final int RHOj3yzz = 19;
    public static final int RHOj3zzz = 20;
    public static final int T1RHOj0 = 21;
    public static final int T2RHOj0 = 22;
    public static final int T3RHOj0 = 23;
    public static final int PHI = 24;
    private static final long serialVersionUID = 1L;
    private final Vector rij = space.makeVector();
    private final Vector rik = space.makeVector();
    private final Vector rkj = space.makeVector();
    private final Vector unitVector = space.makeVector();
    private final Vector sumGiPhi = space.makeVector();
    private final Vector sumGiRhoj0 = space.makeVector();
    private final Vector sumGiRhoj2 = space.makeVector();
    private final Vector sumGiRhoj1x = space.makeVector();
    private final Vector sumGiRhoj1y = space.makeVector();
    private final Vector sumGiRhoj1z = space.makeVector();
    private final Vector sumGiRhoj2xx = space.makeVector();
    private final Vector sumGiRhoj2xy = space.makeVector();
    private final Vector sumGiRhoj2xz = space.makeVector();
    private final Vector sumGiRhoj2yy = space.makeVector();
    private final Vector sumGiRhoj2yz = space.makeVector();
    private final Vector sumGiRhoj2zz = space.makeVector();
    private final Vector sumGiRhoj3xxx = space.makeVector();
    private final Vector sumGiRhoj3xxy = space.makeVector();
    private final Vector sumGiRhoj3xxz = space.makeVector();
    private final Vector sumGiRhoj3xyy = space.makeVector();
    private final Vector sumGiRhoj3xyz = space.makeVector();
    private final Vector sumGiRhoj3xzz = space.makeVector();
    private final Vector sumGiRhoj3yyy = space.makeVector();
    private final Vector sumGiRhoj3yyz = space.makeVector();
    private final Vector sumGiRhoj3yzz = space.makeVector();
    private final Vector sumGiRhoj3zzz = space.makeVector();
    private final Vector sumt1GiRhoj0 = space.makeVector();
    private final Vector sumt2GiRhoj0 = space.makeVector();
    private final Vector sumt3GiRhoj0 = space.makeVector();
    private final Vector rin = space.makeVector();
    private final Vector gjPhi = space.makeVector();
    private final Vector gjRhoj0 = space.makeVector();
    private final Vector gjRhoj1 = space.makeVector();
    private final Vector gjRhoj2 = space.makeVector();
    private final Vector gjRhoj3 = space.makeVector();
    private final Vector gjRhoj1x = space.makeVector();
    private final Vector gjRhoj1y = space.makeVector();
    private final Vector gjRhoj1z = space.makeVector();
    private final Vector gjRhoj2xx = space.makeVector();
    private final Vector gjRhoj2xy = space.makeVector();
    private final Vector gjRhoj2xz = space.makeVector();
    private final Vector gjRhoj2yy = space.makeVector();
    private final Vector gjRhoj2yz = space.makeVector();
    private final Vector gjRhoj2zz = space.makeVector();
    private final Vector gjRhoj3xxx = space.makeVector();
    private final Vector gjRhoj3xxy = space.makeVector();
    private final Vector gjRhoj3xxz = space.makeVector();
    private final Vector gjRhoj3xyy = space.makeVector();
    private final Vector gjRhoj3xyz = space.makeVector();
    private final Vector gjRhoj3xzz = space.makeVector();
    private final Vector gjRhoj3yyy = space.makeVector();
    private final Vector gjRhoj3yyz = space.makeVector();
    private final Vector gjRhoj3yzz = space.makeVector();
    private final Vector gjRhoj3zzz = space.makeVector();
    private final Vector t1GjRhoj0 = space.makeVector();
    private final Vector t2GjRhoj0 = space.makeVector();
    private final Vector t3GjRhoj0 = space.makeVector();
    private final Vector sumGkPhi = space.makeVector();
    private final Vector sumGkRhoj0 = space.makeVector();
    private final Vector sumGkRhoj2 = space.makeVector();
    private final Vector sumGkRhoj1x = space.makeVector();
    private final Vector sumGkRhoj1y = space.makeVector();
    private final Vector sumGkRhoj1z = space.makeVector();
    private final Vector sumGkRhoj2xx = space.makeVector();
    private final Vector sumGkRhoj2xy = space.makeVector();
    private final Vector sumGkRhoj2xz = space.makeVector();
    private final Vector sumGkRhoj2yy = space.makeVector();
    private final Vector sumGkRhoj2yz = space.makeVector();
    private final Vector sumGkRhoj2zz = space.makeVector();
    private final Vector sumGkRhoj3xxx = space.makeVector();
    private final Vector sumGkRhoj3xxy = space.makeVector();
    private final Vector sumGkRhoj3xxz = space.makeVector();
    private final Vector sumGkRhoj3xyy = space.makeVector();
    private final Vector sumGkRhoj3xyz = space.makeVector();
    private final Vector sumGkRhoj3xzz = space.makeVector();
    private final Vector sumGkRhoj3yyy = space.makeVector();
    private final Vector sumGkRhoj3yyz = space.makeVector();
    private final Vector sumGkRhoj3yzz = space.makeVector();
    private final Vector sumGkRhoj3zzz = space.makeVector();
    private final Vector sumt1GkRhoj0 = space.makeVector();
    private final Vector sumt2GkRhoj0 = space.makeVector();
    private final Vector sumt3GkRhoj0 = space.makeVector();
    private final Vector giRij = space.makeVector();
    private final Vector gjRij = space.makeVector();
    private final Vector giRik = space.makeVector();
    private final Vector gjRik = space.makeVector();
    private final Vector giRkj = space.makeVector();
    private final Vector gjRkj = space.makeVector();
    private final Vector giXik = space.makeVector();
    private final Vector gjXik = space.makeVector();
    private final Vector giXkj = space.makeVector();
    private final Vector gjXkj = space.makeVector();
    private final Vector giC = space.makeVector();
    private final Vector gjC = space.makeVector();
    private final Vector giSijk = space.makeVector();
    private final Vector gjSijk = space.makeVector();
    private final Vector giSij = space.makeVector();
    private final Vector gjSij = space.makeVector();
    private final Vector vector100 = space.makeVector();
    private final Vector vector010 = space.makeVector();
    private final Vector vector001 = space.makeVector();
    private final Vector gix = space.makeVector();
    private final Vector giy = space.makeVector();
    private final Vector giz = space.makeVector();
    private final Vector giRhoj0Ref = space.makeVector();
    private final Vector giRhoiRef = space.makeVector();
    private final Vector giEu = space.makeVector();
    private final Vector giFRef = space.makeVector();
    private final Vector giRhoB0 = space.makeVector();
    private final Vector giRhoB2 = space.makeVector();
    private final Vector giRhoSn0 = space.makeVector();
    private final Vector giRhoSn2 = space.makeVector();
    private final Vector giRhoB = space.makeVector();
    private final Vector giRhoSn = space.makeVector();
    private final Vector giFB = space.makeVector();
    private final Vector giFSn = space.makeVector();
    private final Vector giRhoB0Ref = space.makeVector();
    private final Vector giRhoBRef = space.makeVector();
    private final Vector giFBRef = space.makeVector();
    private final Vector giEuB = space.makeVector();
    private final Vector giPhiBB = space.makeVector();
    private final Vector giPhi = space.makeVector();
    private final Vector giRhoj0 = space.makeVector();
    private final Vector giRhoj1 = space.makeVector();
    private final Vector giRhoj2 = space.makeVector();
    private final Vector giRhoj3 = space.makeVector();
    private final Vector giRhoj1x = space.makeVector();
    private final Vector giRhoj1y = space.makeVector();
    private final Vector giRhoj1z = space.makeVector();
    private final Vector giRhoj2xx = space.makeVector();
    private final Vector giRhoj2xy = space.makeVector();
    private final Vector giRhoj2xz = space.makeVector();
    private final Vector giRhoj2yy = space.makeVector();
    private final Vector giRhoj2yz = space.makeVector();
    private final Vector giRhoj2zz = space.makeVector();
    private final Vector giRhoj3xxx = space.makeVector();
    private final Vector giRhoj3xxy = space.makeVector();
    private final Vector giRhoj3xxz = space.makeVector();
    private final Vector giRhoj3xyy = space.makeVector();
    private final Vector giRhoj3xyz = space.makeVector();
    private final Vector giRhoj3xzz = space.makeVector();
    private final Vector giRhoj3yyy = space.makeVector();
    private final Vector giRhoj3yyz = space.makeVector();
    private final Vector giRhoj3yzz = space.makeVector();
    private final Vector giRhoj3zzz = space.makeVector();
    private final Vector t1GiRhoj0 = space.makeVector();
    private final Vector t2GiRhoj0 = space.makeVector();
    private final Vector t3GiRhoj0 = space.makeVector();
    private final Vector gjx = space.makeVector();
    private final Vector gjy = space.makeVector();
    private final Vector gjz = space.makeVector();
    private final Vector ril = space.makeVector();
    private final Vector rlj = space.makeVector();
    private final Vector gkRij = space.makeVector();
    private final Vector gkRik = space.makeVector();
    private final Vector gkRkj = space.makeVector();
    private final Vector gkXik = space.makeVector();
    private final Vector gkXkj = space.makeVector();
    private final Vector gkC = space.makeVector();
    private final Vector gkSijk = space.makeVector();
    private final Vector gkSij = space.makeVector();
    private final Vector gkPhi = space.makeVector();
    private final Vector gkRhoj0 = space.makeVector();
    private final Vector gkRhoj1 = space.makeVector();
    private final Vector gkRhoj2 = space.makeVector();
    private final Vector gkRhoj3 = space.makeVector();
    private final Vector gkRhoj1x = space.makeVector();
    private final Vector gkRhoj1y = space.makeVector();
    private final Vector gkRhoj1z = space.makeVector();
    private final Vector gkRhoj2xx = space.makeVector();
    private final Vector gkRhoj2xy = space.makeVector();
    private final Vector gkRhoj2xz = space.makeVector();
    private final Vector gkRhoj2yy = space.makeVector();
    private final Vector gkRhoj2yz = space.makeVector();
    private final Vector gkRhoj2zz = space.makeVector();
    private final Vector gkRhoj3xxx = space.makeVector();
    private final Vector gkRhoj3xxy = space.makeVector();
    private final Vector gkRhoj3xxz = space.makeVector();
    private final Vector gkRhoj3xyy = space.makeVector();
    private final Vector gkRhoj3xyz = space.makeVector();
    private final Vector gkRhoj3xzz = space.makeVector();
    private final Vector gkRhoj3yyy = space.makeVector();
    private final Vector gkRhoj3yyz = space.makeVector();
    private final Vector gkRhoj3yzz = space.makeVector();
    private final Vector gkRhoj3zzz = space.makeVector();
    private final Vector t1GkRhoj0 = space.makeVector();
    private final Vector t2GkRhoj0 = space.makeVector();
    private final Vector t3GkRhoj0 = space.makeVector();
    private final Vector sumGnPhi = space.makeVector();
    private final Vector gnRhoi0 = space.makeVector();
    private final Vector gnRhoi1sq = space.makeVector();
    private final Vector gnRhoi2sq = space.makeVector();
    private final Vector gnRhoi3sq = space.makeVector();
    private final Vector gntav1 = space.makeVector();
    private final Vector gntav2 = space.makeVector();
    private final Vector gntav3 = space.makeVector();
    private final Vector gnGamma = space.makeVector();
    private final Vector gnRhoi = space.makeVector();
    private final Vector gnF = space.makeVector();
    private final Vector giRhoi0 = space.makeVector();
    private final Vector giRhoi1sq = space.makeVector();
    private final Vector giRhoi2sq = space.makeVector();
    private final Vector giRhoi3sq = space.makeVector();
    private final Vector gitav1 = space.makeVector();
    private final Vector gitav2 = space.makeVector();
    private final Vector gitav3 = space.makeVector();
    private final Vector giGamma = space.makeVector();
    private final Vector giRhoi = space.makeVector();
    private final Vector giF = space.makeVector();
    protected Boundary boundary;
    double jcut = 6.0; //this may not be ideal cutoff for FCC Cu system
    double kcut = jcut * 1.14;
    double[] sum = new double[25];
    private ParameterSetMEAM[] parameters = new ParameterSetMEAM[0];
    private ParameterSetMEAM[] parametersIMC = new ParameterSetMEAM[0];
    private ParameterSetMEAM pi, pj, pk, pl, b, b3sn;
    private ParameterSetMEAM pSn = ParameterSetMEAM.Sn;
    private ParameterSetMEAM pAg = ParameterSetMEAM.Ag;
    private ParameterSetMEAM pCu = ParameterSetMEAM.Cu;
    private Vector[] gnEi = new Vector[0];

    public PotentialMEAM(Space space) {
		super(space);
    }

    public void setParameters(AtomType atomType, ParameterSetMEAM p) {
        int index = atomType.getIndex();
        if(index >= parameters.length) { //15 parameters for each species
			 parameters = (ParameterSetMEAM[])Arrays.resizeArray(parameters, index+1);
		 }
		 parameters[index] = p;
	 }

    public void setParametersIMC(AtomType atomType, ParameterSetMEAM p) {
        //If Cu (Ag) is involved, the reference IMC structure is Cu3Sn (Ag3Sn).
        int index = atomType.getIndex();
        if(index > parametersIMC.length) {
             parametersIMC =
                     (ParameterSetMEAM[])Arrays.resizeArray(parametersIMC, index+1);
		 }
        parametersIMC[index] = p;
    }

	public double getRange() {
		//from MEAMP2
		return kcut;
	}

	public void calcSums(IAtomList atoms) {
		for (int i = 0; i < sum.length; i++) {
    		sum[i] = 0;
		}
        IAtom atom0 = atoms.getAtom(0);
		int indexi = atom0.getType().getIndex(); pi = parameters[indexi];
		for(int j = 1; j < atoms.getAtomCount(); j++) {
            IAtom atomj = atoms.getAtom(j);
			rij.Ev1Mv2(atomj.getPosition(), atom0.getPosition());
			boundary.nearestImage(rij);
            double r = Math.sqrt(rij.squared());
            if (r > jcut) {
                //System.out.println("atom j  "+j+" | rij "+r+", continue.");
                continue;
            }

			int indexj = atomj.getType().getIndex(); pj = parameters[indexj];
            /**To determine amount of screening between atoms i and j
             * by any atom k which may be between them.
			*/
			double Sij = 1.0;
			for(int k = 1; k < atoms.getAtomCount(); k++) {
				if (k == j) continue;
                IAtom atomk = atoms.getAtom(k);
				rik.Ev1Mv2(atomk.getPosition(), atom0.getPosition());
				boundary.nearestImage(rik);
                double ik = Math.sqrt(rik.squared());
                if (ik > r*1.14) continue;
                double v = ((rij.getX(0) * rik.getX(0))
                        +(rij.getX(1) * rik.getX(1))
					        +(rij.getX(2) * rik.getX(2)) ) /(r*ik);
				if (v < -1.0) continue;
				double anglekij = Math.toDegrees(Math.acos(v));
				if (anglekij >= 90) continue;
				rkj.Ev1Mv2(atomk.getPosition(), atomj.getPosition());
				boundary.nearestImage(rkj);
				double kj = Math.sqrt(rkj.squared());
				//System.out.println("	atom k "+k);
				//System.out.println("		rik "+ik+" | rkj "+kj+" | angle "+anglekij);
				//from Baskes (1997)
				double xik = (ik/r)*(ik/r);
				double xjk = (kj/r)*(kj/r);
				double C = ((2.0*(xik + xjk)) - ((xik - xjk)*(xik - xjk))- 1.0)/
					       (1.0 - ((xik - xjk)*(xik - xjk)));
				//System.out.print("		C "+ C);
				if (C < 0) {
					//System.out.println(" | Sijk 1.0 b/c C is negative");
					continue; // negative C forms hyperbola, not ellipse
				}
                int indexk = atomk.getType().getIndex();
                pk = parameters[indexk];

				//Cu-Sn system only

				double Cmin;
				if (pi == pCu & pj == pCu & pk == pCu) Cmin = pCu.Cmin;
				else Cmin = pSn.Cmin;


                //Ag-Sn system only:
				/**
				double Cmin;
				if (pi== pSn & pj == pSn & pk == pSn) Cmin = pSn.Cmin;
				else Cmin = pAg.Cmin;
				*/

				double Sijk;
                if (C <= Cmin) {
                    Sij = 0;
					//System.out.println(" | Sijk 0 ");
					break; // break out of for loop over k atoms
				}
				else if (C >= pi.Cmax) { //Sijk = 1, value of Sij won't change
					//System.out.println(" | Sijk 1.0 ");
        			continue; //continue to next k atom in for loop
        		}
				else {
					double q = ((C - Cmin)/(pi.Cmax - Cmin));
        			Sijk = (1.0 - ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)))
						  *(1.0 - ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)));
        			//System.out.println(" | Sijk " + Sijk);
				}
				Sij *= Sijk;
			} // exit for loop over k atoms

			if (Sij == 0){
			    //System.out.println("Sij is zero for atoms "+indexi+", "+atomj+".");
			    continue;
            } // continue to next j atom in for loop
            //System.out.println("		Sij "+ Sij);

			double rhoj0, rhoj1, rhoj2, rhoj3, x, y, z;

			rhoj0 = pj.rho0 * Math.exp(-pj.b0 * ((r/pj.r0) - 1.0)) * Sij;
			rhoj1 = pj.rho0 * Math.exp(-pj.b1 * ((r/pj.r0) - 1.0)) * Sij;
			rhoj2 = pj.rho0 * Math.exp(-pj.b2 * ((r/pj.r0) - 1.0)) * Sij;
			rhoj3 = pj.rho0 * Math.exp(-pj.b3 * ((r/pj.r0) - 1.0)) * Sij;

			//if(rhoj0==0){System.out.println("rhoj0 is zero.");}

			unitVector.E(rij);
			unitVector.normalize();
			x = unitVector.getX(0);
			y = unitVector.getX(1);
			z = unitVector.getX(2);

            sum[RHOj0] += rhoj0;

			sum[RHOj1x] += rhoj1 * x;
			sum[RHOj1y] += rhoj1 * y;
            sum[RHOj1z] += rhoj1 * z;

			sum[RHOj2xx] += rhoj2 * x * x;
			sum[RHOj2xy] += rhoj2 * x * y;
			sum[RHOj2xz] += rhoj2 * x * z;
			sum[RHOj2yy] += rhoj2 * y * y;
			sum[RHOj2yz] += rhoj2 * y * z;
			sum[RHOj2zz] += rhoj2 * z * z;
			sum[RHOj2] += rhoj2;

			sum[RHOj3xxx] += rhoj3 * x * x * x;
			sum[RHOj3xxy] += rhoj3 * x * x * y;
			sum[RHOj3xxz] += rhoj3 * x * x * z;
			sum[RHOj3xyy] += rhoj3 * x * y * y;
			sum[RHOj3xyz] += rhoj3 * x * y * z;
			sum[RHOj3xzz] += rhoj3 * x * z * z;
			sum[RHOj3yyy] += rhoj3 * y * y * y;
			sum[RHOj3yyz] += rhoj3 * y * y * z;
			sum[RHOj3yzz] += rhoj3 * y * z * z;
			sum[RHOj3zzz] += rhoj3 * z * z * z;

			sum[T1RHOj0] += rhoj0 * pi.t1;
			sum[T2RHOj0] += rhoj0 * pi.t2;
			sum[T3RHOj0] += rhoj0 * pi.t3;

			if (pi.r0 == pj.r0) {
				double a, Eu, rhoj0Ref, rhoiRef, FRef;
				a = pi.a * ((r/pi.r0) - 1.0);
				Eu = - pi.Ec * (1.0 + a) * Math.exp(-a);
				rhoj0Ref = pj.rho0 * Math.exp(-pj.b0 * ((r/pj.r0) - 1.0));
                rhoiRef = pi.Z * rhoj0Ref;
                FRef = pi.A * pi.Ec * (rhoiRef/pi.Z) * Math.log(rhoiRef/pi.Z);
				sum[PHI] += ((2.0/pi.Z) * (Eu - FRef)) * Sij;
			}

			else {
				if (pi == pSn) { // atom i is the Sn atom
					b  = parameters[indexj]; b3sn = parametersIMC[indexj];
				}
				else {
					b  = parameters[indexi]; b3sn = parametersIMC[indexi];
				}
				double a, Eu, rhoB0, rhoB2, rhoSn0, rhoSn2, rhoB, rhoSn, FB, FSn,
					aB, EuB, rhoB0Ref, rhoiBRef, FBRef, phiBB;
				a = b3sn.a * ((r/b3sn.r0) - 1.0);
				Eu = - b3sn.Ec * (1.0 + a) * Math.exp(-a);
				rhoB0 = b.rho0 * Math.exp(-b.b0 * ((r/b.r0) - 1.0));
				rhoB2 = b.rho0 * Math.exp(-b.b2 * ((r/b.r0) - 1.0));
				rhoSn0 = pSn.rho0 * Math.exp(-pSn.b0 * ((r/pSn.r0) - 1.0));
				rhoSn2 = pSn.rho0 * Math.exp(-pSn.b2 * ((r/pSn.r0) - 1.0));
                rhoB = Math.sqrt(
                        ((8.0*rhoB0 + 4.0*rhoSn0)*(8.0*rhoB0 + 4.0*rhoSn0))
					      +((8.0/3.0)*b.t2*(rhoB2 - rhoSn2)*(rhoB2 - rhoSn2)))
					   /(12.0*b.rho0);
				rhoSn = 3.0*rhoB0/pSn.rho0;
				FB  = b.A  * b.Ec  * (rhoB/12.0) * Math.log(rhoB/12.0);
				FSn = pSn.A * pSn.Ec * (rhoSn/4.0) * Math.log(rhoSn/4.0);
				//phiBB
				aB = b.a * ((r/b.r0) - 1.0);
				EuB = - b.Ec * (1.0 + aB) * Math.exp(-aB);
				rhoB0Ref = b.rho0 * Math.exp(-b.b0 * ((r/b.r0) - 1.0));
                rhoiBRef = b.Z * rhoB0Ref;
                FBRef = b.A * b.Ec * (rhoiBRef/b.Z) * Math.log(rhoiBRef/b.Z);
				phiBB = ((2.0/b.Z) * (EuB - FBRef));
				sum[PHI] += ((1.0/3.0)*Eu - 0.25*FB - (1.0/12.0)*FSn - phiBB)*Sij;
			}
		} // exit for loop over j atoms
		//System.out.println("Done");

	} // exit calcSums()

	/** The following methods are not called until after calcSums() is called in
     * energy() and gradient(), so we do not need to call calcSums() in these
     * methods.
	 */

	double rhoi0() {
		return sum[RHOj0]; //
	}

    protected double rhoi1sq() {
    	return    (sum[RHOj1x] * sum[RHOj1x])
	            + (sum[RHOj1y] * sum[RHOj1y])
		   	    + (sum[RHOj1z] * sum[RHOj1z]);
    }

    protected double rhoi2sq() {
    	return (sum[RHOj2xx] * sum[RHOj2xx])
			+ (2.0 * sum[RHOj2xy] * sum[RHOj2xy])
			+ (2.0 * sum[RHOj2xz] * sum[RHOj2xz])
                + (sum[RHOj2yy] * sum[RHOj2yy])
                + (2.0 * sum[RHOj2yz] * sum[RHOj2yz])
                + (sum[RHOj2zz] * sum[RHOj2zz])
                - ((1.0/3.0) * sum[RHOj2] * sum[RHOj2]);
    }

    protected double rhoi3sq() {
    	return     (sum[RHOj3xxx] * sum[RHOj3xxx])
			+(3.0 * sum[RHOj3xxy] * sum[RHOj3xxy])
			+(3.0 * sum[RHOj3xxz] * sum[RHOj3xxz])
			+(3.0 * sum[RHOj3xyy] * sum[RHOj3xyy])
			+(6.0 * sum[RHOj3xyz] * sum[RHOj3xyz])
			+(3.0 * sum[RHOj3xzz] * sum[RHOj3xzz])
			+      (sum[RHOj3yyy] * sum[RHOj3yyy])
			+(3.0 * sum[RHOj3yyz] * sum[RHOj3yyz])
			+(3.0 * sum[RHOj3yzz] * sum[RHOj3yzz])
			+      (sum[RHOj3zzz] * sum[RHOj3zzz]);
    }

    protected double tav1() {
    	return sum[T1RHOj0] / sum[RHOj0];
    }

    protected double tav2() {
    	return sum[T2RHOj0] / sum[RHOj0];
    }

    protected double tav3() {
    	return sum[T3RHOj0] / sum[RHOj0];
    }

    protected double gamma() {
        double rhoi0 = rhoi0(), rhoi1sq = rhoi1sq(), rhoi2sq = rhoi2sq(),
                rhoi3sq = rhoi3sq(), tav1 = tav1(), tav2 = tav2(), tav3 = tav3();
        return ((tav1 * rhoi1sq) + (tav2 * rhoi2sq) + (tav3 * rhoi3sq)) / (rhoi0 * rhoi0);
    }

    protected double rhoi(IAtomList atoms) {
    	double rhoi0 = rhoi0(), gamma = gamma();
		pi = parameters[atoms.getAtom(0).getType().getIndex()];
    	if (pi == pSn) {
    		return (2.0 * rhoi0) / (1.0 + Math.exp(-gamma)); //Sn
    	}
    	return rhoi0 * Math.sqrt(1.0 + gamma); //Cu or Ag
    }

	/* (non-Javadoc)
	 * @see etomica.potential.Potential#energy(etomica.atom.AtomSet)
	 */
	public double energy(IAtomList atoms) {
		calcSums(atoms);
		double rhoi = rhoi(atoms);
		pi = parameters[atoms.getAtom(0).getType().getIndex()];
		double F = pi.A * pi.Ec * (rhoi/pi.Z) * Math.log(rhoi/pi.Z);
		return F + (0.5*sum[PHI]);
	}

	/* (non-Javadoc)
	 * @see etomica.potential.Potential#setBox(etomica.box.Box)
	 */
	public void setBox(Box box) {
		boundary = box.getBoundary();
	}

	public double virial(IAtomList atoms) {
		return calcVirial(atoms, null);
	}

    public Vector[] gradient(IAtomList atoms) {
        calcVirial(atoms, null);
        return gnEi;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        calcVirial(atoms, pressureTensor);
        return gnEi;
    }

    /**
     * This returns the virial (as you'd get from calling virial(atoms)), and
     * fills in the pressureTensor (as you'd get from calling
     * gradient(atoms,pressureTensor)).  It also calculates the gradient but
     * doesn't return it.  virial() and both gradient methods from the
     * PotentialSoft interface call this method and return what the particular
     * method wants.
     */
	private double calcVirial(IAtomList atoms, Tensor pressureTensor) {
        double virial = 0;

		if (atoms.getAtomCount() > gnEi.length) {
			gnEi = new Vector[atoms.getAtomCount()];
			for (int i = 0; i < atoms.getAtomCount(); i++) {
				gnEi[i] = space.makeVector();
			}
		}
		else {
			for (int i = 0; i < atoms.getAtomCount(); i++) {
			    gnEi[i].E(0);
			}
		}

		//check to see if atoms is more than one atom
		if(atoms.getAtomCount()==1){
		    return 0;
		}

        calcSums(atoms);

        if(sum[RHOj0]==0){
            System.out.println("Returning zero for force: atom "+atoms.getAtom(0));
            return 0;
        }

        double rhoi0 = rhoi0(), rhoi1sq = rhoi1sq(), rhoi2sq = rhoi2sq(),
                rhoi3sq = rhoi3sq(), tav1 = tav1(), tav2 = tav2(), tav3 = tav3(),
                gamma = gamma(), rhoi = rhoi(atoms);

        IAtom atom0 = atoms.getAtom(0);
		int indexi = atom0.getType().getIndex(); pi = parameters[indexi];

		sumGiPhi.E(0); sumGiRhoj0.E(0); sumGiRhoj2.E(0);
        sumGiRhoj1x.E(0);
        sumGiRhoj1y.E(0);
        sumGiRhoj1z.E(0);
        sumGiRhoj2xx.E(0);
        sumGiRhoj2xy.E(0);
        sumGiRhoj2xz.E(0);
        sumGiRhoj2yy.E(0); sumGiRhoj2yz.E(0); sumGiRhoj2zz.E(0);
        sumGiRhoj3xxx.E(0);
        sumGiRhoj3xxy.E(0);
        sumGiRhoj3xxz.E(0);
        sumGiRhoj3xyy.E(0);
        sumGiRhoj3xyz.E(0);
        sumGiRhoj3xzz.E(0);
        sumGiRhoj3yyy.E(0);
        sumGiRhoj3yyz.E(0);
        sumGiRhoj3yzz.E(0);
        sumGiRhoj3zzz.E(0);
        sumt1GiRhoj0.E(0);
        sumt2GiRhoj0.E(0);
        sumt3GiRhoj0.E(0);


        for(int n = 1; n < atoms.getAtomCount(); n++) {
            IAtom atomn = atoms.getAtom(n);
            rin.Ev1Mv2(atomn.getPosition(), atom0.getPosition());
            boundary.nearestImage(rin);
            double in = Math.sqrt(rin.squared());
            if (in > kcut) continue; //only consider n that could be j or k to i

            /** We must reset all of the gradients with respect to j, used in
             * calculation of gradient-with-respect-to-n terms, to be the zero
             * vector.  If an atom n is only a k atom to atom i, the gj terms will
             * not be calculated for this n, and, if we don't reset the gj terms
             * to the zero vector as we do below, the gj values for the previous
            * n will be used.
            */

            gjPhi.E(0); gjRhoj0.E(0); gjRhoj2.E(0);
            gjRhoj1x.E(0); gjRhoj1y.E(0); gjRhoj1z.E(0);
            gjRhoj2xx.E(0); gjRhoj2xy.E(0); gjRhoj2xz.E(0);
            gjRhoj2yy.E(0);
            gjRhoj2yz.E(0);
            gjRhoj2zz.E(0);
            gjRhoj3xxx.E(0); gjRhoj3xxy.E(0); gjRhoj3xxz.E(0);
            gjRhoj3xyy.E(0); gjRhoj3xyz.E(0); gjRhoj3xzz.E(0);
            gjRhoj3yyy.E(0); gjRhoj3yyz.E(0); gjRhoj3yzz.E(0); gjRhoj3zzz.E(0);
            t1GjRhoj0.E(0); t2GjRhoj0.E(0); t3GjRhoj0.E(0);

            sumGkPhi.E(0); sumGkRhoj0.E(0); sumGkRhoj2.E(0);
            sumGkRhoj1x.E(0); sumGkRhoj1y.E(0); sumGkRhoj1z.E(0);
            sumGkRhoj2xx.E(0); sumGkRhoj2xy.E(0); sumGkRhoj2xz.E(0);
            sumGkRhoj2yy.E(0);
            sumGkRhoj2yz.E(0);
            sumGkRhoj2zz.E(0);
            sumGkRhoj3xxx.E(0); sumGkRhoj3xxy.E(0); sumGkRhoj3xxz.E(0);
            sumGkRhoj3xyy.E(0); sumGkRhoj3xyz.E(0); sumGkRhoj3xzz.E(0);
            sumGkRhoj3yyy.E(0);
            sumGkRhoj3yyz.E(0);
            sumGkRhoj3yzz.E(0);
            sumGkRhoj3zzz.E(0);
            sumt1GkRhoj0.E(0); sumt2GkRhoj0.E(0); sumt3GkRhoj0.E(0);

            //Here we test to see if n qualifies as a j atom for atom i.
            if (in <= jcut) {
            	rij.E(rin); double ij = in;
    			int indexj = atomn.getType().getIndex(); pj = parameters[indexj];
    			// to calculate Sij, giSij, gjSij
            	double Sij = 1.0; giSij.E(0); gjSij.E(0);
            	for(int k = 1; k < atoms.getAtomCount(); k++) {
                    IAtom atomk = atoms.getAtom(k);
            		if (k == n) continue; // continue to next k atom
            		rik.Ev1Mv2(atomk.getPosition(), atom0.getPosition());
            		boundary.nearestImage(rik);
            		double ik = Math.sqrt(rik.squared());
            		if (ik > ij*1.14) continue; // continue to next k atom
                    double v = ((rij.getX(0) * rik.getX(0))
                            +(rij.getX(1) * rik.getX(1))
					            +(rij.getX(2) * rik.getX(2)) ) /(ij*ik);
            		if (v < -1.0) continue;
            		double anglekij = Math.toDegrees(Math.acos(v));
            		if (anglekij >= 90) continue;
            		rkj.Ev1Mv2(atomk.getPosition(), atomn.getPosition());
            		boundary.nearestImage(rkj);
            		double kj = Math.sqrt(rkj.squared());
            		//from Baskes (1997)
            		double xik = (ik/ij)*(ik/ij);
            		double xkj = (kj/ij)*(kj/ij);
                    double C = ((2.0 * (xik + xkj)) -
                            ((xik - xkj) * (xik - xkj)) - 1.0) /
                            (1.0 - ((xik - xkj)*(xik - xkj)));
            		if (C < 0) continue; // - C does not form ellipse
            		int indexk = atomk.getType().getIndex(); pk = parameters[indexk];

                    //Cu-Sn system only

                    double Cmin;
    				if (pi == pCu & pj == pCu & pk == pCu) {
    					Cmin = pCu.Cmin;
    				}
    				else {
    					Cmin = pSn.Cmin;
    				}


                    //Ag-Sn system only:
    				/**
    				double Cmin;
    				if (pi == pSn & pj == pSn & pk == pSn) Cmin = pSn.Cmin;
    				else Cmin = pAg.Cmin;
    				*/

                    double q = ((C - Cmin)/(pi.Cmax - Cmin));
            		double Sijk;
                    if (C <= Cmin) {
                        Sij = 0;
            			break; //break out of for loop over k atoms
            		}
            		else if (C >= pi.Cmax) { //Sijk = 1, value of Sij won't change
            			continue; // continue to next k atom
            		}
            		else {
            			Sijk = (1.0 - ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)))
    						  *(1.0 - ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)));
            		}

                    giRij.Ea1Tv1(-1.0/ij, rij);
            		gjRij.Ea1Tv1(-1.0, giRij);

                    giRik.Ea1Tv1(-1.0/ik, rik);
            		gjRik.E(0);

                    giRkj.E(0);
            		gjRkj.Ea1Tv1(-1.0/kj, rkj);

                    giXik.Ea1Tv1(-ik/(ij*ij), giRij);
            		giXik.PEa1Tv1(1.0/ij, giRik);
            		giXik.TE(2.0*ik/ij);

                    giXkj.Ea1Tv1(-kj/(ij*ij), giRij);
            		giXkj.PEa1Tv1(1.0/ij, giRkj);
            		giXkj.TE(2.0*kj/ij);

                    gjXik.Ea1Tv1(-ik/(ij*ij), gjRij);
            		gjXik.PEa1Tv1(1.0/ij, gjRik);
            		gjXik.TE(2.0*ik/ij);

                    gjXkj.Ea1Tv1(-kj/(ij*ij), gjRij);
            		gjXkj.PEa1Tv1(1.0/ij, gjRkj);
            		gjXkj.TE(2.0*kj/ij);

                    giC.Ea1Tv1( 1.0 + (xik - xkj)*(C - 1.0), giXik);
            		giC.PEa1Tv1(1.0 - (xik - xkj)*(C + 1.0), giXkj);
            		giC.TE(2.0/(1.0 - ((xik - xkj)*(xik - xkj))));

                    gjC.Ea1Tv1( 1.0 + (xik - xkj)*(C - 1.0), gjXik);
            		gjC.PEa1Tv1(1.0 - (xik - xkj)*(C + 1.0), gjXkj);
            		gjC.TE(2.0/(1.0 - ((xik - xkj)*(xik - xkj))));

                    giSijk.Ea1Tv1(8
                            * (1.0 - ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)))
				            * ((1.0 - q)*(1.0 - q)*(1.0 - q))
							* (1.0 / (pi.Cmax - Cmin)), giC);

                    gjSijk.Ea1Tv1(8
                            * (1.0 - ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)))
		                    * ((1.0 - q)*(1.0 - q)*(1.0 - q))
							* (1.0 / (pi.Cmax - Cmin)), gjC);

                    /** The Sij value used to calculate gradSij is that for
                     * previous k's, or, for the first k considered, 1.0.  The
                     * same goes for the giSij/gjSij values, except these are
                     * initialized as the zero vector...
            		 */

                    giSij.TE(Sijk);
            		giSij.PEa1Tv1(Sij, giSijk);

                    gjSij.TE(Sijk);
            		gjSij.PEa1Tv1(Sij, gjSijk);

                    Sij *= Sijk;
            	} // exit loop over k for n = j

                if (Sij != 0) {

                    double rhoj0, rhoj1, rhoj2, rhoj3, x, y, z, phi;

                    rhoj0 = pj.rho0 * Math.exp(-pj.b0 * ((ij/pj.r0) - 1.0)) * Sij;
	            	rhoj1 = pj.rho0 * Math.exp(-pj.b1 * ((ij/pj.r0) - 1.0)) * Sij;
	            	rhoj2 = pj.rho0 * Math.exp(-pj.b2 * ((ij/pj.r0) - 1.0)) * Sij;
	            	rhoj3 = pj.rho0 * Math.exp(-pj.b3 * ((ij/pj.r0) - 1.0)) * Sij;

                    unitVector.E(rij);
	            	unitVector.normalize();
	            	x = unitVector.getX(0);
	            	y = unitVector.getX(1);
	            	z = unitVector.getX(2);

                    vector100.setX(0,1.0);
	        		vector010.setX(1,1.0);
	        		vector001.setX(2,1.0);

                    gix.Ea1Tv1(x/(ij*ij), rij);
	        		gix.PEa1Tv1(-1.0/ij, vector100);

                    giy.Ea1Tv1(y/(ij*ij), rij);
	        		giy.PEa1Tv1(-1.0/ij, vector010);

                    giz.Ea1Tv1(z/(ij*ij), rij);
	        		giz.PEa1Tv1(-1.0/ij, vector001);

                    giRij.Ea1Tv1(-1.0/ij, rij);

                    if (pi == pj) {
	        			double a, Eu, rhoj0Ref, rhoiRef, FRef;
	    				a = pi.a * ((ij/pi.r0) - 1.0);
	    				Eu = - pi.Ec * (1.0 + a) * Math.exp(-a);
	    				rhoj0Ref = pj.rho0 * Math.exp(-pj.b0 * ((ij/pj.r0) - 1.0));
                        rhoiRef = pi.Z * rhoj0Ref;
                        FRef = pi.A * pi.Ec * (rhoiRef/pi.Z) * Math.log(rhoiRef/pi.Z);
	    				phi = ((2.0/pi.Z) * (Eu - FRef)) * Sij;

                        giRhoj0Ref.Ea1Tv1(-rhoj0Ref*pj.b0/pj.r0, giRij);
	    			    giRhoiRef.Ea1Tv1(pi.Z, giRhoj0Ref);
	    			    giFRef.Ea1Tv1( (pi.A*pi.Ec/pi.Z)
	    			        	      *(1.0 + Math.log(rhoiRef/pi.Z)), giRhoiRef);
                        giEu.Ea1Tv1((pi.Ec * pi.a * pi.a / pi.r0)
                                        *((ij/pi.r0) - 1.0)
                                        * (Math.exp(-pi.a * ((ij / pi.r0) - 1.0))),
                                giRij);
                        giPhi.E(giEu);
	    			    giPhi.ME(giFRef);
	    			    giPhi.TE(2.0 * Sij /pi.Z);
	    			    giPhi.PEa1Tv1(phi/Sij, giSij);
	    			    sumGiPhi.PE(giPhi);

                        gjPhi.Ea1Tv1(-1.0, giEu);
		        		gjPhi.PE(giFRef);
		        		gjPhi.TE(2.0 * Sij /pi.Z);
		        		gjPhi.PEa1Tv1(phi/Sij, gjSij);
	    			} else {
	    				if (pi == pSn) { // atom i is the Sn atom
	    					b  = parameters[indexj]; b3sn = parametersIMC[indexj];
	    				}
	    				else {
	    					b  = parameters[indexi]; b3sn = parametersIMC[indexi];
	    				}
	    				double a, Eu, rhoB0, rhoB2, rhoSn0, rhoSn2, rhoB, rhoSn,
							FB, FSn, aB, EuB, rhoB0Ref, rhoiBRef, FBRef, phiBB;
	    				a = b3sn.a * ((ij/b3sn.r0) - 1.0);
	    				Eu = - b3sn.Ec * (1.0 + a) * Math.exp(-a);
	    				rhoB0 = b.rho0 * Math.exp(-b.b0 * ((ij/b.r0) - 1.0));
	    				rhoB2 = b.rho0 * Math.exp(-b.b2 * ((ij/b.r0) - 1.0));
	    				rhoSn0 = pSn.rho0 * Math.exp(-pSn.b0 * ((ij/pSn.r0) - 1.0));
	    				rhoSn2 = pSn.rho0 * Math.exp(-pSn.b2 * ((ij/pSn.r0) - 1.0));
                        rhoB = Math.sqrt(
                                ((8.0*rhoB0 + 4.0*rhoSn0)*(8.0*rhoB0 + 4.0*rhoSn0))
	    					    +((8.0/3.0)*b.t2*(rhoB2 - rhoSn2)*(rhoB2 - rhoSn2)))
	    					   /(12.0*b.rho0);
	    				rhoSn = 3.0*rhoB0/pSn.rho0;
	    				FB  = b.A  * b.Ec  * (rhoB/12.0) * Math.log(rhoB/12.0);
	    				FSn = pSn.A * pSn.Ec * (rhoSn/4.0) * Math.log(rhoSn/4.0);
	    				//phiBB
	    				aB = b.a * ((ij/b.r0) - 1.0);
	    				EuB = - b.Ec * (1.0 + aB) * Math.exp(-aB);
	    				rhoB0Ref = b.rho0 * Math.exp(-b.b0 * ((ij/b.r0) - 1.0));
                        rhoiBRef = b.Z * rhoB0Ref;
                        FBRef = b.A * b.Ec * (rhoiBRef/b.Z) * Math.log(rhoiBRef/b.Z);
	    				phiBB = ((2.0/b.Z) * (EuB - FBRef));
	    				phi = ((1.0/3.0)*Eu - 0.25*FB - (1.0/12.0)*FSn - phiBB)*Sij;

                        //giPhi
	    				giRhoB0.Ea1Tv1(-rhoB0*b.b0/b.r0, giRij);
	    				giRhoB0.Ea1Tv1(-rhoB2*b.b2/b.r0, giRij);
	    				giRhoSn0.Ea1Tv1(-rhoSn0*pSn.b0/pSn.r0, giRij);
	    				giRhoSn0.Ea1Tv1(-rhoSn2*pSn.b2/pSn.r0, giRij);

                        giRhoB.Ea1Tv1(16.0*(8.0*rhoB0 + 4.0*rhoSn0), giRhoB0);
	    				giRhoB.PEa1Tv1(16.0*(8.0*rhoB0 + 4.0*rhoSn0), giRhoSn0);
	    				giRhoB.PEa1Tv1( (16.0/3.0)*b.t2*(rhoB2 - rhoSn2), giRhoB2);
	    				giRhoB.PEa1Tv1(-(16.0/3.0)*b.t2*(rhoB2 - rhoSn2), giRhoSn2);
	    				giRhoB.TE(1.0/(2.0*144.0*b.rho0*b.rho0*rhoB));

                        giRhoSn.Ea1Tv1(3.0/pSn.rho0, giRhoB0);
	        			giFB.Ea1Tv1(  (b.A * b.Ec / b.Z)
	        					     *(1.0 + Math.log(rhoB/12.0)), giRhoB);
	        			giFSn.Ea1Tv1( (pSn.A * pSn.Ec / pSn.Z)
	        					     *(1.0 + Math.log(rhoSn/4.0)), giRhoSn);
                        giEu.Ea1Tv1((b3sn.Ec * b3sn.a * b3sn.a / b3sn.r0)
                                        *((ij/b3sn.r0) - 1.0)
                                        * (Math.exp(-b3sn.a * ((ij / b3sn.r0) - 1.0))),
                                giRij);
	        			//giPhiBB
	        			giRhoB0Ref.Ea1Tv1(-rhoB0Ref*b.b0/b.r0, giRij);
	    			    giRhoBRef.Ea1Tv1(b.Z, giRhoB0Ref);
	    			    giFBRef.Ea1Tv1( (b.A*b.Ec/b.Z)
	    			        	       *(1.0 + Math.log(rhoiBRef/pi.Z)), giRhoiRef);
	    			    giEuB.Ea1Tv1( (b.Ec*b.a*b.a/b.r0) * ((ij/b.r0) - 1.0)
	    			        	   	   *(Math.exp(-b.a*((ij/b.r0)-1.0))), giRij);
	    			    giPhiBB.E(giEuB);
	    			    giPhiBB.ME(giFBRef);
	    			    giPhiBB.TE(2.0/b.Z);

                        giPhi.Ea1Tv1(1.0/3.0, giEu);
	        			giPhi.PEa1Tv1(-0.25, giFB);
	        			giPhi.PEa1Tv1(-1.0/12.0, giFSn);
	        			giPhi.PEa1Tv1(-1.0, giPhiBB);
	        			giPhi.TE(Sij);
	    			    giPhi.PEa1Tv1(phi/Sij, giSij);
	    			    sumGiPhi.PE(giPhi);

                        gjPhi.Ea1Tv1(-1.0/3.0, giEu);
	        			gjPhi.PEa1Tv1(0.25, giFB);
	        			gjPhi.PEa1Tv1(1.0/12.0, giFSn);
	        			gjPhi.PEa1Tv1(1.0, giPhiBB);
		        		gjPhi.TE(Sij);
		        		gjPhi.PEa1Tv1(phi/Sij, gjSij);
	        		}

	        		giRhoj0.Ea1Tv1(rhoj0/Sij, giSij);
	        		giRhoj0.PEa1Tv1(-rhoj0*pj.b0/(pj.r0), giRij);
	        		sumGiRhoj0.PE(giRhoj0);

                    giRhoj1.Ea1Tv1(rhoj1/Sij, giSij);
	        		giRhoj1.PEa1Tv1(-rhoj1*pj.b1/(pj.r0), giRij);

                    giRhoj1x.Ea1Tv1(rhoj1, gix);
	        		giRhoj1x.PEa1Tv1(x, giRhoj1);
	        		sumGiRhoj1x.PE(giRhoj1x);

                    giRhoj1y.Ea1Tv1(rhoj1, giy);
	        		giRhoj1y.PEa1Tv1(y, giRhoj1);
	        		sumGiRhoj1y.PE(giRhoj1y);

                    giRhoj1z.Ea1Tv1(rhoj1, giz);
	        		giRhoj1z.PEa1Tv1(z, giRhoj1);
	        		sumGiRhoj1z.PE(giRhoj1z);

                    giRhoj2.Ea1Tv1(rhoj2/Sij, giSij);
	        		giRhoj2.PEa1Tv1(-rhoj2*pj.b2/(pj.r0), giRij);
	        		sumGiRhoj2.PE(giRhoj2);

                    giRhoj2xx.Ea1Tv1(2.0*rhoj2*x, gix);
	        		giRhoj2xx.PEa1Tv1(x*x, giRhoj2);
	        		sumGiRhoj2xx.PE(giRhoj2xx);

                    giRhoj2xy.Ea1Tv1(rhoj2*x, giy);
	        		giRhoj2xy.PEa1Tv1(rhoj2*y, gix);
	        		giRhoj2xy.PEa1Tv1(x*y, giRhoj2);
	        		sumGiRhoj2xy.PE(giRhoj2xy);

                    giRhoj2xz.Ea1Tv1(rhoj2*x, giz);
	        		giRhoj2xz.PEa1Tv1(rhoj2*z, gix);
	        		giRhoj2xz.PEa1Tv1(x*z, giRhoj2);
	        		sumGiRhoj2xz.PE(giRhoj2xz);

	        		giRhoj2yy.Ea1Tv1(2.0*rhoj2*y, giy);
	        		giRhoj2yy.PEa1Tv1(y*y, giRhoj2);
	        		sumGiRhoj2yy.PE(giRhoj2yy);

                    giRhoj2yz.Ea1Tv1(rhoj2*y, giz);
	        		giRhoj2yz.PEa1Tv1(rhoj2*z, giy);
	        		giRhoj2yz.PEa1Tv1(y*z, giRhoj2);
	        		sumGiRhoj2yz.PE(giRhoj2yz);

                    giRhoj2zz.Ea1Tv1(2.0*rhoj2*z, giz);
	        		giRhoj2zz.PEa1Tv1(z*z, giRhoj2);
	        		sumGiRhoj2zz.PE(giRhoj2zz);

                    giRhoj3.Ea1Tv1(rhoj3/Sij, giSij);
	        		giRhoj3.PEa1Tv1(-rhoj3*pj.b3/(pj.r0), giRij);

                    giRhoj3xxx.Ea1Tv1(3.0*rhoj3*x*x, gix);
	        		giRhoj3xxx.PEa1Tv1(x*x*x, giRhoj3);
	        		sumGiRhoj3xxx.PE(giRhoj3xxx);

                    giRhoj3xxy.Ea1Tv1(rhoj3*x*x, giy);
	        		giRhoj3xxy.PEa1Tv1(2.0*rhoj3*x*y, gix);
	        		giRhoj3xxy.PEa1Tv1(x*x*y, giRhoj3);
	        		sumGiRhoj3xxy.PE(giRhoj3xxy);

	        		giRhoj3xxz.Ea1Tv1(rhoj3*x*x, giz);
	        		giRhoj3xxz.PEa1Tv1(2.0*rhoj3*x*z, gix);
	        		giRhoj3xxz.PEa1Tv1(x*x*z, giRhoj3);
	        		sumGiRhoj3xxz.PE(giRhoj3xxz);

                    giRhoj3xyy.Ea1Tv1(2.0*rhoj3*x*y, giy);
	        		giRhoj3xyy.PEa1Tv1(rhoj3*y*y, gix);
	        		giRhoj3xyy.PEa1Tv1(x*y*y, giRhoj3);
	        		sumGiRhoj3xyy.PE(giRhoj3xyy);

                    giRhoj3xyz.Ea1Tv1(rhoj3*x*y, giz);
	        		giRhoj3xyz.PEa1Tv1(rhoj3*x*z, giy);
	        		giRhoj3xyz.PEa1Tv1(rhoj3*y*z, gix);
	        		giRhoj3xyz.PEa1Tv1(x*y*z, giRhoj3);
	        		sumGiRhoj3xyz.PE(giRhoj3xyz);

                    giRhoj3xzz.Ea1Tv1(2.0*rhoj3*x*z, giz);
	        		giRhoj3xzz.PEa1Tv1(rhoj3*z*z, gix);
	        		giRhoj3xzz.PEa1Tv1(x*z*z, giRhoj3);
	        		sumGiRhoj3xzz.PE(giRhoj3xzz);

	        		giRhoj3yyy.Ea1Tv1(3.0*rhoj3*y*y, giy);
	        		giRhoj3yyy.PEa1Tv1(y*y*y, giRhoj3);
	        		sumGiRhoj3yyy.PE(giRhoj3yyy);

                    giRhoj3yyz.Ea1Tv1(rhoj3*y*y, giz);
	        		giRhoj3yyz.PEa1Tv1(2.0*rhoj3*y*z, giy);
	        		giRhoj3yyz.PEa1Tv1(y*y*z, giRhoj3);
	        		sumGiRhoj3yyz.PE(giRhoj3yyz);

                    giRhoj3yzz.Ea1Tv1(2.0*rhoj3*y*z, giz);
	        		giRhoj3yzz.PEa1Tv1(rhoj3*z*z, giy);
	        		giRhoj3yzz.PEa1Tv1(y*z*z, giRhoj3);
	        		sumGiRhoj3yzz.PE(giRhoj3yzz);

                    giRhoj3zzz.Ea1Tv1(3.0*rhoj3*z*z, giz);
	        		giRhoj3zzz.PEa1Tv1(z*z*z, giRhoj3);
	        		sumGiRhoj3zzz.PE(giRhoj3zzz);

                    t1GiRhoj0.Ea1Tv1(pi.t1, giRhoj0);
	        		sumt1GiRhoj0.PE(t1GiRhoj0);

                    t2GiRhoj0.Ea1Tv1(pi.t2, giRhoj0);
	        		sumt2GiRhoj0.PE(t2GiRhoj0);

                    t3GiRhoj0.Ea1Tv1(pi.t3, giRhoj0);
	        		sumt3GiRhoj0.PE(t3GiRhoj0);

                    //Pair-wise terms required to calculate gnEi
	        		gjx.Ea1Tv1(-1.0, gix);
	        		gjy.Ea1Tv1(-1.0, giy);
	        		gjz.Ea1Tv1(-1.0, giz);
	        		gjRij.Ea1Tv1(-1.0, giRij);

                    gjRhoj0.Ea1Tv1(rhoj0/Sij, gjSij);
	        		gjRhoj0.PEa1Tv1(-rhoj0*pj.b0/(pj.r0), gjRij);

                    gjRhoj1.Ea1Tv1(rhoj1/Sij, gjSij);
	        		gjRhoj1.PEa1Tv1(-rhoj1*pj.b1/(pj.r0), gjRij);

                    gjRhoj1x.Ea1Tv1(rhoj1, gjx);
	        		gjRhoj1x.PEa1Tv1(x, gjRhoj1);

                    gjRhoj1y.Ea1Tv1(rhoj1, gjy);
	        		gjRhoj1y.PEa1Tv1(y, gjRhoj1);

                    gjRhoj1z.Ea1Tv1(rhoj1, gjz);
	        		gjRhoj1z.PEa1Tv1(z, gjRhoj1);

                    gjRhoj2.Ea1Tv1(rhoj2/Sij, gjSij);
	        		gjRhoj2.PEa1Tv1(-rhoj2*pj.b2/(pj.r0), gjRij);

                    gjRhoj2xx.Ea1Tv1(2.0*rhoj2*x, gjx);
	        		gjRhoj2xx.PEa1Tv1(x*x, gjRhoj2);

                    gjRhoj2xy.Ea1Tv1(rhoj2*x, gjy);
	        		gjRhoj2xy.PEa1Tv1(rhoj2*y, gjx);
	        		gjRhoj2xy.PEa1Tv1(x*y, gjRhoj2);

                    gjRhoj2xz.Ea1Tv1(rhoj2*x, gjz);
	        		gjRhoj2xz.PEa1Tv1(rhoj2*z, gjx);
	        		gjRhoj2xz.PEa1Tv1(x*z, gjRhoj2);

                    gjRhoj2yy.Ea1Tv1(2.0*rhoj2*y, gjy);
	        		gjRhoj2yy.PEa1Tv1(y*y, gjRhoj2);

                    gjRhoj2yz.Ea1Tv1(rhoj2*y, gjz);
	        		gjRhoj2yz.PEa1Tv1(rhoj2*z, gjy);
	        		gjRhoj2yz.PEa1Tv1(y*z, gjRhoj2);

                    gjRhoj2zz.Ea1Tv1(2.0*rhoj2*z, gjz);
	        		gjRhoj2zz.PEa1Tv1(z*z, gjRhoj2);

                    gjRhoj3.Ea1Tv1(rhoj3/Sij, gjSij);
	        		gjRhoj3.PEa1Tv1(-rhoj3*pj.b3/(pj.r0), gjRij);

                    gjRhoj3xxx.Ea1Tv1(3.0*rhoj3*x*x, gjx);
	        		gjRhoj3xxx.PEa1Tv1(x*x*x, gjRhoj3);

                    gjRhoj3xxy.Ea1Tv1(rhoj3*x*x, gjy);
	        		gjRhoj3xxy.PEa1Tv1(2.0*rhoj3*x*y, gjx);
	        		gjRhoj3xxy.PEa1Tv1(x*x*y, gjRhoj3);

                    gjRhoj3xxz.Ea1Tv1(rhoj3*x*x, gjz);
	        		gjRhoj3xxz.PEa1Tv1(2.0*rhoj3*x*z, gjx);
	        		gjRhoj3xxz.PEa1Tv1(x*x*z, gjRhoj3);

                    gjRhoj3xyy.Ea1Tv1(2.0*rhoj3*x*y, gjy);
	        		gjRhoj3xyy.PEa1Tv1(rhoj3*y*y, gjx);
	        		gjRhoj3xyy.PEa1Tv1(x*y*y, gjRhoj3);

                    gjRhoj3xyz.Ea1Tv1(rhoj3*x*y, gjz);
	        		gjRhoj3xyz.PEa1Tv1(rhoj3*x*z, gjy);
	        		gjRhoj3xyz.PEa1Tv1(rhoj3*y*z, gjx);
	        		gjRhoj3xyz.PEa1Tv1(x*y*z, gjRhoj3);

                    gjRhoj3xzz.Ea1Tv1(2.0*rhoj3*x*z, gjz);
	        		gjRhoj3xzz.PEa1Tv1(rhoj3*z*z, gjx);
	        		gjRhoj3xzz.PEa1Tv1(x*z*z, gjRhoj3);

                    gjRhoj3yyy.Ea1Tv1(3.0*rhoj3*y*y, gjy);
	        		gjRhoj3yyy.PEa1Tv1(y*y*y, gjRhoj3);

                    gjRhoj3yyz.Ea1Tv1(rhoj3*y*y, gjz);
	        		gjRhoj3yyz.PEa1Tv1(2.0*rhoj3*y*z, gjy);
	        		gjRhoj3yyz.PEa1Tv1(y*y*z, gjRhoj3);

                    gjRhoj3yzz.Ea1Tv1(2.0*rhoj3*y*z, gjz);
	        		gjRhoj3yzz.PEa1Tv1(rhoj3*z*z, gjy);
	        		gjRhoj3yzz.PEa1Tv1(y*z*z, gjRhoj3);

                    gjRhoj3zzz.Ea1Tv1(3.0*rhoj3*z*z, gjz);
	        		gjRhoj3zzz.PEa1Tv1(z*z*z, gjRhoj3);

                    t1GjRhoj0.Ea1Tv1(pi.t1, gjRhoj0);

                    t2GjRhoj0.Ea1Tv1(pi.t2, gjRhoj0);

                    t3GjRhoj0.Ea1Tv1(pi.t3, gjRhoj0);
            	} // exit if statement with condition that Sij != 0
            } // exit if statement with condition that n is a j atom of i

            /** To consider n as a k atom, we must loop through neighbors j of i
             * again.
             */

            for(int j = 1; j < atoms.getAtomCount(); j++) {
            	//The k atom, n, must not be treated as one of the other j atoms.
            	if (j == n) continue; // continue to next j atom
                IAtom atomj = atoms.getAtom(j);
        		rij.Ev1Mv2(atomj.getPosition(), atom0.getPosition());
        		boundary.nearestImage(rij);
        		double ij = Math.sqrt(rij.squared());
        		if (ij > jcut) continue; // continue to next j atom
        		rik.E(rin); double ik = in;
	        	if (ik > ij*1.14) continue; // n won't impact this i-j interaction
                double v = ((rij.getX(0) * rik.getX(0))
                        +(rij.getX(1) * rik.getX(1))
				            +(rij.getX(2) * rik.getX(2)) ) /(ij*ik);
	        	if (v < -1.0) continue;
	        	double anglekij = Math.toDegrees(Math.acos(v));
	        	if (anglekij >= 90) continue;
	        	rkj.Ev1Mv2(atomn.getPosition(), atomj.getPosition());
	        	boundary.nearestImage(rkj);
	        	double kj = Math.sqrt(rkj.squared());
	        	//from Baskes (1997)
	        	double xik = (ik/ij)*(ik/ij);
	        	double xkj = (kj/ij)*(kj/ij);
	        	double C = ( (2.0*(xik + xkj)) - ((xik - xkj)*(xik - xkj))- 1.0 )
							/ (1.0 - ((xik - xkj)*(xik - xkj)));
	        	if (C < 0) continue;
	        	int indexj = atomj.getType().getIndex(); pj = parameters[indexj];
	        	int indexk = atomn.getType().getIndex(); pk = parameters[indexk];

                //Cu-Sn system only

                double Cmin;
				if (pi == pCu & pj == pCu & pk == pCu) Cmin = pCu.Cmin;
				else Cmin = pSn.Cmin;


                //Ag-Sn system only:
				/**
				double Cmin;
				if (pi == pSn & pj == pSn & pk == pSn) Cmin = pSn.Cmin;
				else Cmin = pAg.Cmin;
			    */

                double q = ((C - Cmin)/(pi.Cmax - Cmin));
	        	double Sijk;
                if (C <= Cmin) {
                    continue;
                    /** If Sijk equals 0 or 1, any gradient of Sijk,
                     * including that with respect to k, will also be 0.  We
                     * should continue on to the next j atom, whose gkSijk may be
                     * nonzero.
                     */
	        	} else if (C >= pi.Cmax) {
                    continue; //see above reasoning
	        	}
	        	else {
        			Sijk = (1.0 - ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)))
						  *(1.0 - ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)));
	        	}

                // To calculate Sij. l is k != n.
	        	// We can start Sij out with the value for Sijk (k = n).
	        	double Sij = Sijk;
	        	for(int l = 1; l < atoms.getAtomCount(); l++) {
	        		if (l == j || l == n) continue; //already have Sijk for n = k
                    IAtom atoml = atoms.getAtom(l);
	    			ril.Ev1Mv2(atoml.getPosition(), atom0.getPosition());
	    			boundary.nearestImage(ril);
	    			double il = Math.sqrt(ril.squared());
	    			if (il > ij*1.14) continue;
                    double w = ((rij.getX(0) * ril.getX(0))
                            +(rij.getX(1) * ril.getX(1))
					            +(rij.getX(2) * ril.getX(2)) ) /(ij*il);
	    			if (w < -1.0) continue;
	    			double anglelij = Math.toDegrees(Math.acos(w));
	    			if (anglelij >= 90) continue;
	    			rlj.Ev1Mv2(atoml.getPosition(), atomj.getPosition());
	    			boundary.nearestImage(rlj);
	    			double lj = Math.sqrt(rlj.squared());
	    			//from Baskes (1997)
	    			double xil = (il/ij)*(il/ij);
	    			double xjl = (lj/ij)*(lj/ij);
                    double c = ((2.0 * (xil + xjl)) -
                            ( (xil - xjl)*(xil - xjl) ) - 1.0 ) /
	    					   ( 1.0 - ( (xil - xjl)*(xil - xjl) ) );
	    			if (c < 0) continue;
	    			int indexl = atoml.getType().getIndex(); pl = parameters[indexl];

                    //Cu-Sn system only

                    double cmin;
    				if (pi == pCu & pj == pCu & pl == pCu) cmin = pCu.Cmin;
    				else cmin = pSn.Cmin;


                    //Ag-Sn system only:
    				/**
    				double cmin;
    				if (pi == pSn & pj == pSn & pk == pSn) cmin = pSn.Cmin;
    				else cmin = pAg.Cmin;
    				*/

                    double Sijl;
                    if (c <= cmin) {
                        Sij = 0;
	    				break; // break out of loop over k atoms for j!=n atom
	    			}
	    			else if (c >= pi.Cmax) { //Sijl = 1, value of Sij won't change
	            		continue; // continue to next k!=n atom
	            	}
	    			else {
	    				double m = ((c - cmin)/(pi.Cmax - cmin));
	        			Sijl = (1.0 - ((1.0 - m)*(1.0 - m)*(1.0 - m)*(1.0 - m)))
							  *(1.0 - ((1.0 - m)*(1.0 - m)*(1.0 - m)*(1.0 - m)));
	    			}
	    			Sij *= Sijl;
                } // exit loop over k!=n for j!=n

	        	// An l (k!=n) atom may have made Sij = 0.
	        	if (Sij == 0) continue; //continue to next j!=n

                double rhoj0, rhoj1, rhoj2, rhoj3, x, y, z, phi;

                rhoj0 = pj.rho0 * Math.exp(-pj.b0 * ((ij/pj.r0) - 1.0)) * Sij;
	        	rhoj1 = pj.rho0 * Math.exp(-pj.b1 * ((ij/pj.r0) - 1.0)) * Sij;
	            rhoj2 = pj.rho0 * Math.exp(-pj.b2 * ((ij/pj.r0) - 1.0)) * Sij;
	        	rhoj3 = pj.rho0 * Math.exp(-pj.b3 * ((ij/pj.r0) - 1.0)) * Sij;

                unitVector.E(rij);
	            unitVector.normalize();
	            x = unitVector.getX(0);
	            y = unitVector.getX(1);
	            z = unitVector.getX(2);

                if (pi == pj) {
        			double a, EuRef, rhoj0Ref, rhoiRef, FRef;
    				a = pi.a * ((ij/pi.r0) - 1.0);
    				EuRef = - pi.Ec * (1.0 + a) * Math.exp(-a);
    				rhoj0Ref = pj.rho0 * Math.exp(-pj.b0 * ((ij/pj.r0) - 1.0));
                    rhoiRef = pi.Z * rhoj0Ref;
                    FRef = pi.A * pi.Ec * (rhoiRef/pi.Z) * Math.log(rhoiRef/pi.Z);
    				phi = ((2.0/pi.Z) * (EuRef - FRef)) * Sij;
    			} else {
    				if (pi == pSn) { // atom i is the Sn atom
    					b  = parameters[indexj]; b3sn = parametersIMC[indexj];
    				}
    				else {
    					b  = parameters[indexi]; b3sn = parametersIMC[indexi];
    				}
    				double a, Eu, rhoB0, rhoB2, rhoSn0, rhoSn2, rhoB, rhoSn, FB,
						FSn, aB, EuB, rhoB0Ref, rhoiBRef, FBRef, phiBB;
    				a = b3sn.a * ((ij/b3sn.r0) - 1.0);
    				Eu = - b3sn.Ec * (1.0 + a) * Math.exp(-a);
    				rhoB0 = b.rho0 * Math.exp(-b.b0 * ((ij/b.r0) - 1.0));
    				rhoB2 = b.rho0 * Math.exp(-b.b2 * ((ij/b.r0) - 1.0));
    				rhoSn0 = pSn.rho0 * Math.exp(-pSn.b0 * ((ij/pSn.r0) - 1.0));
    				rhoSn2 = pSn.rho0 * Math.exp(-pSn.b2 * ((ij/pSn.r0) - 1.0));
                    rhoB = Math.sqrt(
                            ((8.0*rhoB0 + 4.0*rhoSn0)*(8.0*rhoB0 + 4.0*rhoSn0))
    					     +((8.0/3.0)*b.t2*(rhoB2 - rhoSn2)*(rhoB2 - rhoSn2)))
    					   /(12.0*b.rho0);
    				rhoSn = 3.0*rhoB0/pSn.rho0;
    				FB  = b.A  * b.Ec  * (rhoB/12.0) * Math.log(rhoB/12.0);
    				FSn = pSn.A * pSn.Ec * (rhoSn/4.0) * Math.log(rhoSn/4.0);
    				//phiBB
    				aB = b.a * ((ij/b.r0) - 1.0);
    				EuB = - b.Ec * (1.0 + aB) * Math.exp(-aB);
    				rhoB0Ref = b.rho0 * Math.exp(-b.b0 * ((ij/b.r0) - 1.0));
                    rhoiBRef = b.Z * rhoB0Ref;
                    FBRef = b.A * b.Ec * (rhoiBRef/b.Z) * Math.log(rhoiBRef/b.Z);
    				phiBB = ((2.0/b.Z) * (EuB - FBRef));
    				phi = ((1.0/3.0)*Eu - 0.25*FB - (1.0/12.0)*FSn - phiBB) * Sij;
        		}

                gkRij.E(0);
                gkRik.Ea1Tv1(1.0 / ik, rik);
                gkRkj.Ea1Tv1(1.0 / kj, rkj);

	        	gkXik.Ea1Tv1(-ik/(ij*ij), gkRij);
	        	gkXik.PEa1Tv1(1.0/ij, gkRik);
	        	gkXik.TE(2.0*ik/ij);

                gkXkj.Ea1Tv1(-kj/(ij*ij), gkRij);
	        	gkXkj.PEa1Tv1(1.0/ij, gkRkj);
	        	gkXkj.TE(2.0*kj/ij);

                gkC.Ea1Tv1( 1.0 + (xik - xkj)*(C - 1.0), gkXik);
		    	gkC.PEa1Tv1(1.0 - (xik - xkj)*(C + 1.0), gkXkj);
		    	gkC.TE( 2.0 / ( 1.0 - ((xik - xkj)*(xik - xkj)) ));

                gkSijk.Ea1Tv1(8 * (1.0 - ((1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)))
	                            * ((1.0 - q)*(1.0 - q)*(1.0 - q))
						        * (1.0 / (pi.Cmax - Cmin)), gkC);

                //We only consider one k atom - the k atom that is n
		    	gkSij.Ea1Tv1(Sij/Sijk, gkSijk);

                gkPhi.Ea1Tv1(phi/Sij, gkSij); sumGkPhi.PE(gkPhi);

                gkRhoj0.Ea1Tv1(rhoj0/Sij, gkSij); sumGkRhoj0.PE(gkRhoj0);
		    	gkRhoj1.Ea1Tv1(rhoj1/Sij, gkSij);
		    	gkRhoj2.Ea1Tv1(rhoj2/Sij, gkSij); sumGkRhoj2.PE(gkRhoj2);
		    	gkRhoj3.Ea1Tv1(rhoj3/Sij, gkSij);

                gkRhoj1x.Ea1Tv1(x, gkRhoj1); sumGkRhoj1x.PE(gkRhoj1x);
		    	gkRhoj1y.Ea1Tv1(y, gkRhoj1); sumGkRhoj1y.PE(gkRhoj1y);
		    	gkRhoj1z.Ea1Tv1(z, gkRhoj1); sumGkRhoj1z.PE(gkRhoj1z);

                gkRhoj2xx.Ea1Tv1(x*x, gkRhoj2); sumGkRhoj2xx.PE(gkRhoj2xx);
		    	gkRhoj2xy.Ea1Tv1(x*y, gkRhoj2); sumGkRhoj2xy.PE(gkRhoj2xy);
		    	gkRhoj2xz.Ea1Tv1(x*z, gkRhoj2); sumGkRhoj2xz.PE(gkRhoj2xz);
		    	gkRhoj2yy.Ea1Tv1(y*y, gkRhoj2); sumGkRhoj2yy.PE(gkRhoj2yy);
		    	gkRhoj2yz.Ea1Tv1(y*z, gkRhoj2); sumGkRhoj2yz.PE(gkRhoj2yz);
		    	gkRhoj2zz.Ea1Tv1(z*z, gkRhoj2); sumGkRhoj2zz.PE(gkRhoj2zz);

                gkRhoj3xxx.Ea1Tv1(x*x*x, gkRhoj3); sumGkRhoj3xxx.PE(gkRhoj3xxx);
		    	gkRhoj3xxy.Ea1Tv1(x*x*y, gkRhoj3); sumGkRhoj3xxy.PE(gkRhoj3xxy);
		    	gkRhoj3xxz.Ea1Tv1(x*x*z, gkRhoj3); sumGkRhoj3xxz.PE(gkRhoj3xxz);
		    	gkRhoj3xyy.Ea1Tv1(x*y*y, gkRhoj3); sumGkRhoj3xyy.PE(gkRhoj3xyy);
		    	gkRhoj3xyz.Ea1Tv1(x*y*z, gkRhoj3); sumGkRhoj3xyz.PE(gkRhoj3xyz);
		        gkRhoj3xzz.Ea1Tv1(x*z*z, gkRhoj3); sumGkRhoj3xzz.PE(gkRhoj3xzz);
		    	gkRhoj3yyy.Ea1Tv1(y*y*y, gkRhoj3); sumGkRhoj3yyy.PE(gkRhoj3yyy);
		    	gkRhoj3yyz.Ea1Tv1(y*y*z, gkRhoj3); sumGkRhoj3yyz.PE(gkRhoj3yyz);
		    	gkRhoj3yzz.Ea1Tv1(y*z*z, gkRhoj3); sumGkRhoj3yzz.PE(gkRhoj3yzz);
		    	gkRhoj3zzz.Ea1Tv1(z*z*z, gkRhoj3); sumGkRhoj3zzz.PE(gkRhoj3zzz);

                t1GkRhoj0.Ea1Tv1(pi.t1, gkRhoj0); sumt1GkRhoj0.PE(t1GkRhoj0);
		    	t2GkRhoj0.Ea1Tv1(pi.t2, gkRhoj0); sumt2GkRhoj0.PE(t2GkRhoj0);
		    	t3GkRhoj0.Ea1Tv1(pi.t3, gkRhoj0); sumt3GkRhoj0.PE(t3GkRhoj0);
	        } //exit loop over j!=n atoms, with n as a k atom

            //multi-body terms, n as k is included

            sumGnPhi.E(sumGkPhi);
	    	sumGnPhi.PE(gjPhi);

            gnRhoi0.E(sumGkRhoj0);
	    	gnRhoi0.PE(gjRhoj0);

            gnRhoi1sq.Ea1Tv1( sum[RHOj1x], sumGkRhoj1x);
	    	gnRhoi1sq.PEa1Tv1(sum[RHOj1x], gjRhoj1x);
			gnRhoi1sq.PEa1Tv1(sum[RHOj1y], sumGkRhoj1y);
			gnRhoi1sq.PEa1Tv1(sum[RHOj1y], gjRhoj1y);
			gnRhoi1sq.PEa1Tv1(sum[RHOj1z], sumGkRhoj1z);
			gnRhoi1sq.PEa1Tv1(sum[RHOj1z], gjRhoj1z);
			gnRhoi1sq.TE(2.0);

            gnRhoi2sq.Ea1Tv1( 2.0 * sum[RHOj2xx], sumGkRhoj2xx);
			gnRhoi2sq.PEa1Tv1(2.0 * sum[RHOj2xx], gjRhoj2xx);
			gnRhoi2sq.PEa1Tv1(4.0 * sum[RHOj2xy], sumGkRhoj2xy);
			gnRhoi2sq.PEa1Tv1(4.0 * sum[RHOj2xy], gjRhoj2xy);
			gnRhoi2sq.PEa1Tv1(4.0 * sum[RHOj2xz], sumGkRhoj2xz);
			gnRhoi2sq.PEa1Tv1(4.0 * sum[RHOj2xz], gjRhoj2xz);
			gnRhoi2sq.PEa1Tv1(2.0 * sum[RHOj2yy], sumGkRhoj2yy);
			gnRhoi2sq.PEa1Tv1(2.0 * sum[RHOj2yy], gjRhoj2yy);
			gnRhoi2sq.PEa1Tv1(4.0 * sum[RHOj2yz], sumGkRhoj2yz);
			gnRhoi2sq.PEa1Tv1(4.0 * sum[RHOj2yz], gjRhoj2yz);
			gnRhoi2sq.PEa1Tv1(2.0 * sum[RHOj2zz], sumGkRhoj2zz);
			gnRhoi2sq.PEa1Tv1(2.0 * sum[RHOj2zz], gjRhoj2zz);
			gnRhoi2sq.PEa1Tv1(-(2.0/3.0) * sum[RHOj2], sumGkRhoj2);
			gnRhoi2sq.PEa1Tv1(-(2.0/3.0) * sum[RHOj2], gjRhoj2);


            gnRhoi3sq.Ea1Tv1( sum[RHOj3xxx], sumGkRhoj3xxx);
			gnRhoi3sq.PEa1Tv1(sum[RHOj3xxx], gjRhoj3xxx);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xxy], sumGkRhoj3xxy);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xxy], gjRhoj3xxy);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xxz], sumGkRhoj3xxz);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xxz], gjRhoj3xxz);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xyy], sumGkRhoj3xyy);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xyy], gjRhoj3xyy);
			gnRhoi3sq.PEa1Tv1(6.0 * sum[RHOj3xyz], sumGkRhoj3xyz);
			gnRhoi3sq.PEa1Tv1(6.0 * sum[RHOj3xyz], gjRhoj3xyz);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xzz], sumGkRhoj3xzz);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xzz], gjRhoj3xzz);
			gnRhoi3sq.PEa1Tv1(sum[RHOj3yyy], sumGkRhoj3yyy);
			gnRhoi3sq.PEa1Tv1(sum[RHOj3yyy], gjRhoj3yyy);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3yyz], sumGkRhoj3yyz);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3yyz], gjRhoj3yyz);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3yzz], sumGkRhoj3yzz);
			gnRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3yzz], gjRhoj3yzz);
			gnRhoi3sq.PEa1Tv1(sum[RHOj3zzz], sumGkRhoj3zzz);
			gnRhoi3sq.PEa1Tv1(sum[RHOj3zzz], gjRhoj3zzz);
			gnRhoi3sq.TE(2.0);

            gntav1.Ea1Tv1(-sum[T1RHOj0]/(rhoi0*rhoi0), gnRhoi0);
			gntav1.PEa1Tv1(1.0/rhoi0, sumt1GkRhoj0);
			gntav1.PEa1Tv1(1.0/rhoi0, t1GjRhoj0);

            gntav2.Ea1Tv1(-sum[T2RHOj0]/(rhoi0*rhoi0), gnRhoi0);
			gntav2.PEa1Tv1(1.0/rhoi0, sumt2GkRhoj0);
			gntav2.PEa1Tv1(1.0/rhoi0, t2GjRhoj0);

            gntav3.Ea1Tv1(-sum[T3RHOj0]/(rhoi0*rhoi0), gnRhoi0);
			gntav3.PEa1Tv1(1.0/rhoi0, sumt3GkRhoj0);
			gntav3.PEa1Tv1(1.0/rhoi0, t3GjRhoj0);

            gnGamma.Ea1Tv1((-2.0/rhoi0)
				*((tav1*rhoi1sq) + (tav2*rhoi2sq) + (tav3*rhoi3sq)), gnRhoi0);
			gnGamma.PEa1Tv1(tav1, gnRhoi1sq);
			gnGamma.PEa1Tv1(tav2, gnRhoi2sq);
			gnGamma.PEa1Tv1(tav3, gnRhoi3sq);
			gnGamma.PEa1Tv1(rhoi1sq, gntav1);
			gnGamma.PEa1Tv1(rhoi2sq, gntav2);
			gnGamma.PEa1Tv1(rhoi3sq, gntav3);
			gnGamma.TE(1.0/(rhoi0*rhoi0));

            if (pi == pSn){
				gnRhoi.Ea1Tv1(rhoi0*
						Math.exp(-gamma)/(1.0 + Math.exp(-gamma)), gnGamma);
				gnRhoi.PE(gnRhoi0);
				gnRhoi.TE(2.0/(1.0+Math.exp(-gamma)));
			}
			else {
				gnRhoi.Ea1Tv1(rhoi0*0.5*Math.sqrt(1.0/(1.0 + gamma)), gnGamma);
				gnRhoi.PEa1Tv1(Math.sqrt(1.0+gamma), gnRhoi0);
			}

            gnF.Ea1Tv1( (pi.A*pi.Ec/pi.Z)*
					(1.0 + Math.log(rhoi/pi.Z)), gnRhoi);

            //System.out.println("gradF is " + gradF);
			//System.exit(0);

            gnEi[n].E(gnF);
            gnEi[n].PEa1Tv1(0.5, sumGnPhi);

            // calculate contribution to pressure tensor and virial
            if (pressureTensor != null) {
                // pressure tensor might be null if we're actually looking for
                // just the virial or just the gradient
                pressureTensor.PEv1v2(gnEi[n], rin);
            }

            virial += gnEi[n].dot(rin);
        } //exit loop over atom n

        giRhoi0.E(sumGiRhoj0);

        giRhoi1sq.Ea1Tv1( sum[RHOj1x], sumGiRhoj1x);
		giRhoi1sq.PEa1Tv1(sum[RHOj1y], sumGiRhoj1y);
		giRhoi1sq.PEa1Tv1(sum[RHOj1z], sumGiRhoj1z);
		giRhoi1sq.TE(2.0);

        giRhoi2sq.Ea1Tv1( 2.0 * sum[RHOj2xx], sumGiRhoj2xx);
		giRhoi2sq.PEa1Tv1(4.0 * sum[RHOj2xy], sumGiRhoj2xy);
		giRhoi2sq.PEa1Tv1(4.0 * sum[RHOj2xz], sumGiRhoj2xz);
		giRhoi2sq.PEa1Tv1(2.0 * sum[RHOj2yy], sumGiRhoj2yy);
		giRhoi2sq.PEa1Tv1(4.0 * sum[RHOj2yz], sumGiRhoj2yz);
		giRhoi2sq.PEa1Tv1(2.0 * sum[RHOj2zz], sumGiRhoj2zz);
		giRhoi2sq.PEa1Tv1(-(2.0/3.0) * sum[RHOj2], sumGiRhoj2);

        giRhoi3sq.Ea1Tv1(sum[RHOj3xxx], sumGiRhoj3xxx);
		giRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xxy], sumGiRhoj3xxy);
		giRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xxz], sumGiRhoj3xxz);
		giRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xyy], sumGiRhoj3xyy);
		giRhoi3sq.PEa1Tv1(6.0 * sum[RHOj3xyz], sumGiRhoj3xyz);
		giRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3xzz], sumGiRhoj3xzz);
		giRhoi3sq.PEa1Tv1(sum[RHOj3yyy], sumGiRhoj3yyy);
		giRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3yyz], sumGiRhoj3yyz);
		giRhoi3sq.PEa1Tv1(3.0 * sum[RHOj3yzz], sumGiRhoj3yzz);
		giRhoi3sq.PEa1Tv1(sum[RHOj3zzz], sumGiRhoj3zzz);
		giRhoi3sq.TE(2.0);

        gitav1.Ea1Tv1(-sum[T1RHOj0]/(rhoi0*rhoi0), sumGiRhoj0);
		gitav1.PEa1Tv1(1.0/rhoi0, sumt1GiRhoj0);

        gitav2.Ea1Tv1(-sum[T2RHOj0]/(rhoi0*rhoi0), sumGiRhoj0);
		gitav2.PEa1Tv1(1.0/rhoi0, sumt2GiRhoj0);

        gitav3.Ea1Tv1(-sum[T3RHOj0]/(rhoi0*rhoi0), sumGiRhoj0);
		gitav3.PEa1Tv1(1.0/rhoi0, sumt3GiRhoj0);

        giGamma.Ea1Tv1((-2.0/rhoi0)
				*((tav1*rhoi1sq) + (tav2*rhoi2sq) + (tav3*rhoi3sq)), giRhoi0);
		giGamma.PEa1Tv1(tav1, giRhoi1sq);
		giGamma.PEa1Tv1(tav2, giRhoi2sq);
		giGamma.PEa1Tv1(tav3, giRhoi3sq);
		giGamma.PEa1Tv1(rhoi1sq, gitav1);
		giGamma.PEa1Tv1(rhoi2sq, gitav2);
		giGamma.PEa1Tv1(rhoi3sq, gitav3);
		giGamma.TE(1.0/(rhoi0*rhoi0));

        if (pi == pSn) {
			giRhoi.Ea1Tv1(rhoi0*Math.exp(-gamma)/(1.0 + Math.exp(-gamma)), giGamma);
			giRhoi.PE(giRhoi0);
			giRhoi.TE(2.0/(1.0+Math.exp(-gamma)));
		}

		else {
			giRhoi.Ea1Tv1(rhoi0*0.5*Math.sqrt(1.0/(1.0 + gamma)), giGamma);
			giRhoi.PEa1Tv1(Math.sqrt(1.0+gamma), giRhoi0);
		}

        giF.Ea1Tv1( (pi.A*pi.Ec/pi.Z)*
				(1.0 + Math.log(rhoi/pi.Z)), giRhoi);

        //System.out.println("gradF is " + gradF);
		//System.exit(0);

        gnEi[0].E(giF);
		gnEi[0].PEa1Tv1(0.5, sumGiPhi);

        if(gnEi[0].isNaN()){
            throw new RuntimeException("atom " + indexi + "    " + atoms.getAtom(indexi));
        }

        return virial;
	}
}
