/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;
import etomica.space.Vector;

public class PotentialCalculationLJSP implements PotentialCalculation {
		
	public PotentialCalculationLJSP(Space space, Box box, CoordinateDefinition coordinateDefinition, double temperature, double[] elasticParams){
		sum = new double[65];
        this.space = space;
        this.temperature = temperature;
        rij = space.makeVector();
        Rij = space.makeVector();
        drj = space.makeVector();
        dri = space.makeVector();
        drij = space.makeVector();

		gx1  = elasticParams[2];  gy1  = elasticParams[3];  gy4  = elasticParams[4];
		gx11 = elasticParams[5];  gy11 = elasticParams[6];  gx44 = elasticParams[7];
		gy44 = elasticParams[8];  gx12 = elasticParams[9];  gz12 = elasticParams[10];

		hx11 = gx11 + gx1*gx1 + 1.0;
		hy11 = gy11 + gy1*gy1;
		hx12 = gx12 + gx1*gy1; //gx2=gy1
		hz12 = gz12 + gy1*gy1; //gz1=gy1 , gz2=gy1
		hx44 = gx44; //gx4=0
		hy44 = gy44 + gy4*gy4 + 1.0/4.0;

		this.box = box;
        this.coordinateDefinition = coordinateDefinition;
        volume = coordinateDefinition.getBox().getBoundary().volume();
        nMol = coordinateDefinition.getBox().getLeafList().size();
	}
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
		Potential2SoftSpherical potentialSoft = (Potential2SoftSpherical)potential;
        Vector[] g = potentialSoft.gradient(atoms);// gradient do nearestImage() !!!
        int nNbrAtoms = atoms.size();
        Vector ri = atoms.get(0).getPosition();
        Vector Ri = coordinateDefinition.getLatticePosition(atoms.get(0));
        dri.Ev1Mv2(ri, Ri);
		box.getBoundary().nearestImage(dri);
		Vector rj;
        
        for (int j=1;j<nNbrAtoms;j++){//START from "1" NOT "0" because we need j != i
        	rj = atoms.get(j).getPosition();
        	Vector Rj = coordinateDefinition.getLatticePosition(atoms.get(j));
        	Rij.Ev1Mv2(Ri , Rj);
        	box.getBoundary().nearestImage(Rij);
			drj.Ev1Mv2(rj , Rj);
			box.getBoundary().nearestImage(drj);
			drij.Ev1Mv2(dri , drj);
			rij.Ev1Pv2(Rij, drij);
			double rij2 = rij.squared();
			double Rij2 = Rij.squared();

			Vector drij_1 = space.makeVector();
			drij_1.setX(0, (gx1 -1.0)*drij.getX(0));
			drij_1.setX(1, gy1 *drij.getX(1));
			drij_1.setX(2, gy1 *drij.getX(2));

			Vector drij_2 = space.makeVector();
			drij_2.setX(0, gy1 *drij.getX(0));
			drij_2.setX(1, (gx1 -1.0)*drij.getX(1));
			drij_2.setX(2, gy1 *drij.getX(2));

			Vector drij_3 = space.makeVector();
			drij_3.setX(0, gy1 *drij.getX(0));
			drij_3.setX(1, gy1 *drij.getX(1));
			drij_3.setX(2, (gx1 -1)*drij.getX(2));

			Vector drij_4 = space.makeVector();
			drij_4.setX(1, drij.getX(2));
			drij_4.setX(2, drij.getX(1));
			drij_4.TE(gy4 - 0.5);

			Vector drij_5 = space.makeVector();
			drij_5.setX(0, drij.getX(2));
			drij_5.setX(2, drij.getX(0));
			drij_5.TE(gy4 - 0.5);

			Vector drij_6 = space.makeVector();
			drij_6.setX(0, drij.getX(1));
			drij_6.setX(1, drij.getX(0));
			drij_6.TE(gy4 - 0.5);

			Vector drij_11 = space.makeVector();
			drij_11.setX(0, hx11*drij.getX(0));
			drij_11.setX(1, hy11*drij.getX(1));
			drij_11.setX(2, hy11*drij.getX(2));

			Vector drij_22 = space.makeVector();
			drij_22.setX(0, hy11*drij.getX(0));
			drij_22.setX(1, hx11*drij.getX(1));
			drij_22.setX(2, hy11*drij.getX(2));

			Vector drij_33 = space.makeVector();
			drij_33.setX(0, hy11*drij.getX(0));
			drij_33.setX(1, hy11*drij.getX(1));
			drij_33.setX(2, hx11*drij.getX(2));

			Vector drij_12 = space.makeVector();
			drij_12.setX(0, hx12*drij.getX(0));
			drij_12.setX(1, hx12*drij.getX(1));
			drij_12.setX(2, hz12*drij.getX(2));

			Vector drij_13 = space.makeVector();
			drij_13.setX(0, hx12*drij.getX(0));
			drij_13.setX(1, hz12*drij.getX(1));
			drij_13.setX(2, hx12*drij.getX(2));

			Vector drij_23 = space.makeVector();
			drij_23.setX(0, hz12*drij.getX(0));
			drij_23.setX(1, hx12*drij.getX(1));
			drij_23.setX(2, hx12*drij.getX(2));

			Vector drij_44 = space.makeVector(); //yz
			drij_44.setX(0, hx44*drij.getX(0));
			drij_44.setX(1, hy44*drij.getX(1));
			drij_44.setX(2, hy44*drij.getX(2));

			Vector drij_55 = space.makeVector();//xz
			drij_55.setX(0, hy44*drij.getX(0));
			drij_55.setX(1, hx44*drij.getX(1));
			drij_55.setX(2, hy44*drij.getX(2));

			Vector drij_66 = space.makeVector();//xy
			drij_66.setX(0, hy44*drij.getX(0));
			drij_66.setX(1, hy44*drij.getX(1));
			drij_66.setX(2, hx44*drij.getX(2));


			Vector rijx1 = space.makeVector();
			rijx1.setX(0, rij.getX(0));
			Vector rijx2 = space.makeVector();
			rijx2.setX(1, rij.getX(0));
			Vector rijx3 = space.makeVector();
			rijx3.setX(2, rij.getX(0));

			Vector rijy1 = space.makeVector();
			rijy1.setX(0, rij.getX(1));
			Vector rijy2 = space.makeVector();
			rijy2.setX(1, rij.getX(1));
			Vector rijy3 = space.makeVector();
			rijy3.setX(2, rij.getX(1));

			Vector rijz1 = space.makeVector();
			rijz1.setX(0, rij.getX(2));
			Vector rijz2 = space.makeVector();
			rijz2.setX(1, rij.getX(2));
			Vector rijz3 = space.makeVector();
			rijz3.setX(2, rij.getX(2));

			Vector rijy3z2 = space.makeVector();
			rijy3z2.Ev1Pv2(rijy3, rijz2);
			Vector rijx3z1 = space.makeVector();
			rijx3z1.Ev1Pv2(rijx3, rijz1);
			Vector rijx2y1 = space.makeVector();
			rijx2y1.Ev1Pv2(rijx2, rijy1);

			//b_mn
			Vector drij_1G = space.makeVector();
			drij_1G.setX(0, gx1 *drij.getX(0));
			drij_1G.setX(1, gy1 *drij.getX(1));
			drij_1G.setX(2, gy1 *drij.getX(2));

			Vector drij_2G = space.makeVector();
			drij_2G.setX(0, gy1 *drij.getX(0));
			drij_2G.setX(1, gx1 *drij.getX(1));
			drij_2G.setX(2, gy1 *drij.getX(2));

			Vector drij_3G = space.makeVector();
			drij_3G.setX(0, gy1 *drij.getX(0));
			drij_3G.setX(1, gy1 *drij.getX(1));
			drij_3G.setX(2, gx1 *drij.getX(2));


			double dW  = potentialSoft.du(rij2);
	        double d2W = potentialSoft.d2u(rij2);

			sum[0] += g[j].dot(rij );  //Fr
	        sum[1] += g[j].dot(drij); //Fdr
			sum[2] += r1Phir2( rij, rij , dW, d2W); //rPhir
			sum[3] += r1Phir2(drij, drij, dW, d2W); //drPhidr
			sum[4] += r1Phir2( rij, drij, dW, d2W); //rPhidr

			sum[5] += g[j].dot(drij_1); //Fdr1
			sum[6] += g[j].dot(drij_2); //Fdr2
			sum[7] += g[j].dot(drij_3); //Fdr3
			sum[8]  += g[j].dot(drij_4); //Fdr4
			sum[9]  += g[j].dot(drij_5); //Fdr5
			sum[10] += g[j].dot(drij_6); //Fdr6

			sum[11] += g[j].dot(drij_11); //Fdr11
			sum[12] += g[j].dot(drij_22); //Fdr22
			sum[13] += g[j].dot(drij_33); //Fdr33

			sum[14] += g[j].dot(drij_12); //Fdr12
			sum[15] += g[j].dot(drij_13); //Fdr13
			sum[16] += g[j].dot(drij_23); //Fdr23

			sum[17] += g[j].dot(drij_44); //Fdr44
			sum[18] += g[j].dot(drij_55); //Fdr55
			sum[19] += g[j].dot(drij_66); //Fdr66

			sum[20]  += g[j].getX(0)* rij.getX(0);//Fxrx
			sum[21]  += g[j].getX(1)* rij.getX(1);//Fyry
			sum[22]  += g[j].getX(2)* rij.getX(2);//Fyrz

			sum[23]  += g[j].getX(0)* rij.getX(1); //Fxry:
			sum[24]  += g[j].getX(0)* rij.getX(2); //Fxrz:
			sum[25]  += g[j].getX(1)* rij.getX(2); //Fyrz:

			sum[26] += rxPhi2ry(0, 0, rij.getX(0), rij.getX(0) ,  dW, d2W); //x_Phixx_x
			sum[27] += rxPhi2ry(1, 1, rij.getX(1), rij.getX(1) ,  dW, d2W); //y_Phiyy_y
			sum[28] += rxPhi2ry(2, 2, rij.getX(2), rij.getX(2) ,  dW, d2W); //z_Phizz_z

			sum[29] += rxPhi2ry(0, 1, rij.getX(0), rij.getX(1) ,  dW, d2W); //x_Phixy_y
			sum[30] += rxPhi2ry(0, 2, rij.getX(0), rij.getX(2) ,  dW, d2W); //x_Phixz_z
			sum[31] += rxPhi2ry(1, 2, rij.getX(1), rij.getX(2) ,  dW, d2W); //y_Phiyz_z

			sum[32] += rxPhi2ry(1, 0, rij.getX(0), rij.getX(1) ,  dW, d2W); //x_Phiyx_y
			sum[33] += rxPhi2ry(2, 0, rij.getX(0), rij.getX(2) ,  dW, d2W); //x_Phizx_z
			sum[34] += rxPhi2ry(2, 1, rij.getX(1), rij.getX(2) ,  dW, d2W); //y_Phizy_z

			sum[35] += r1Phir2(drij_1, drij_1, dW, d2W); //dr1_Phi_dr1
			sum[36] += r1Phir2(drij_2, drij_2, dW, d2W); //dr2_Phi_dr2
			sum[37] += r1Phir2(drij_3, drij_3, dW, d2W); //dr3_Phi_dr3

			sum[38] += r1Phir2(drij_1, drij_2, dW, d2W); //dr1_Phi_dr2
			sum[39] += r1Phir2(drij_1, drij_3, dW, d2W); //dr1_Phi_dr3
			sum[40] += r1Phir2(drij_2, drij_3, dW, d2W); //dr2_Phi_dr3

			sum[41] += r1Phir2(drij_4, drij_4, dW, d2W); //dr4_Phi_dr4
			sum[42] += r1Phir2(drij_5, drij_5, dW, d2W); //dr5_Phi_dr5
			sum[43] += r1Phir2(drij_6, drij_6, dW, d2W); //dr6_Phi_dr6

			sum[44] += r1Phir2(rijx1, drij_1, dW, d2W); //rx1_Phi_dr1
			sum[45] += r1Phir2(rijy2, drij_2, dW, d2W); //ry2_Phi_dr2
			sum[46] += r1Phir2(rijz3, drij_3, dW, d2W); //rz3_Phi_dr3

			sum[47] += r1Phir2(rijx1, drij_2, dW, d2W); //rx1_Phi_dr2
			sum[48] += r1Phir2(rijx1, drij_3, dW, d2W); //rx1_Phi_dr3
			sum[49] += r1Phir2(rijy2, drij_3, dW, d2W); //ry2_Phi_dr3

			sum[50] += r1Phir2(rijy2, drij_1, dW, d2W); //ry2_Phi_dr1
			sum[51] += r1Phir2(rijz3, drij_1, dW, d2W); //rz3_Phi_dr1
			sum[52] += r1Phir2(rijz3, drij_2, dW, d2W); //rz3_Phi_dr2

			sum[53] += r1Phir2(rijy3z2 , drij_4, dW, d2W); //ry3z2_Phi_dr4
			sum[54] += r1Phir2(rijx3z1 , drij_5, dW, d2W); //rx3z1_Phi_dr5
			sum[55] += r1Phir2(rijx2y1 , drij_6, dW, d2W); //rx2y1_Phi_dr6

			//b_mn
			sum[56] += g[j].dot(drij_1G); //Fdr1G
			sum[57] += g[j].dot(drij_2G); //Fdr2G
			sum[58] += g[j].dot(drij_3G); //Fdr3G

			sum[59] += r1Phir2(drij, drij_1, dW, d2W); //drPhidr1
			sum[60] += r1Phir2(drij, drij_2, dW, d2W); //drPhidr2
			sum[61] += r1Phir2(drij, drij_3, dW, d2W); //drPhidr3

			sum[62] += r1Phir2(rijx1, drij, dW, d2W); //rx1_Phi_dr
			sum[63] += r1Phir2(rijy2, drij, dW, d2W); //ry2_Phi_dr
			sum[64] += r1Phir2(rijz3, drij, dW, d2W); //rz3_Phi_dr
		}//j atoms loop
	}

	protected double r1Phir2(Vector T1ij, Vector T2ij, double dW, double d2W){
		double rij2 = rij.squared();
		double rDr = dW/rij2 * T1ij.dot(T2ij) + (d2W-dW)/rij2/rij2 * T1ij.dot(rij)*(T2ij.dot(rij));
		return rDr;
	}

	protected double rxPhi2ry(int a, int b, double t1, double t2 , double dW, double d2W){
		double rij2 = rij.squared();
		double rDr;
		rDr = (d2W-dW)/rij2/rij2 *t1*t2* rij.getX(a)*rij.getX(b);
		return rDr;
	}

	/**
	 * Sets the virial sum to zero, typically to begin a new virial-sum calculation.
	 * @return this instance, so the method can be called in-line as the instance is
	 * passed to the PotentialMaster.
	 */
	public void reset() {
		for(int i=0;i<sum.length;i++){
			sum[i] = 0.0;
		}
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double[] getSum() {
//		System.out.println(sum[2] +"  ,   "+ sum[3]+">>>>>"+(sum[2]+sum[3]));
		return sum;
	}
	
	private double[] sum;
	protected double volume , temperature;
	protected double gx1, gy1, gy4, gx11, gy11, gx12, gz12, gx44, gy44;
	protected double hx11, hy11, hx12, hz12, hx44, hy44;
	protected int nMol;
    protected final Vector rij;
    protected final Vector Rij;
    protected final Vector drj;
    protected final Vector dri ;
    protected final Vector drij;
    protected final Space space;
    protected final Box box;
    protected final CoordinateDefinition coordinateDefinition;
 }
