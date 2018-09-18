/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.molecule.IMoleculeList;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.PotentialMaster;

/**
 *  Determine the second derivative of the atomic/ molecular potential energy w.r.t. to
 *   its generalized coordinates, u, where u is defined as the relative deviation of the 
 *   atom/ molecule from its nominal position 
 * 
 *  The class use <CoordinateDefinition> class SetToU method to put the atom/ molecule in
 *   space and the calculate the change in potential energy. 
 *   
 *  Output is the second derivative of the energy, phi: d2_phi/du_i^2, or d2_phi/(du_i.du_j)
 *  For Hessian Matrix:   d2_phi/(du_i.du_j) = d2_phi/(du_j.du_i)
 * 
 *  NOTE: fixedDeltaU in the field is the interval of change. It has to be CAREFULLY picked!
 *        Usually, start from a larger value (~0.1, e.g.) and the algorithm will reduce 
 *        the value. One way to check is to compare d2f/(duiduj) and d2f/(dujdui), of course
 *        the smaller the value the better it is.
 * 
 * @author taitan
 *
 */
public class CalcNumerical2ndDerivative{
	
	public CalcNumerical2ndDerivative(Box box, PotentialMaster potentialMaster,
                                      CoordinateDefinition coordinateDefinition){
		this.coordinateDefinition = coordinateDefinition;
		meterPotential = new MeterPotentialEnergy(potentialMaster, box);
		
		a = new double[ntab][ntab];
	}
 	
	public double f(double[] newU) {
	
		for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, newU);
		}
		
		return meterPotential.getDataAsScalar();
	}

	public double dphi_du(int[] d, double u[]) {
		
		/*
		 * d[i] = the i-th elements of the derivative
		 * u[i] = the generalized coordinate/ relative displacement from 
		 * 			nominal position 
		 */
        deltaU = fixedDeltaU;

        a[0][0] = computeA(d, u, deltaU);
     
		double err = big;
        double dfdu = Double.NaN;
		
		for(int i=1; i<ntab; i++){
			deltaU = deltaU /con;
            a[0][i] = computeA(d, u, deltaU);
            
			fac = con2;
			
			for(int j=1; j<=i; j++){
				a[j][i] = (a[j-1][i]*fac - a[j-1][i-1])/(fac-1);
				fac = con2*fac;
				errt = Math.max(Math.abs(a[j][i]-a[j-1][i]), Math.abs(a[j][i]-a[j-1][i-1]));
				
				if (errt <= err){
					err = errt;
					dfdu = a[j][i];
				}
			}
			
			if (Math.abs(a[i][i]-a[i-1][i-1]) >= safe*err){
				break;
			}
		}
		
		return dfdu;
	}
	
	public double d2phi_du2(int[] d, double u[]) {
		
		/*
		 * d[i] = the i-th elements of the derivative
		 * u[i] = the generalized coordinate/ relative displacement from 
		 * 			nominal position 
		 */
        deltaU = fixedDeltaU;
        a[0][0] = computeA(d, u, deltaU);
        
		double err = big;
        double d2fdu2 = Double.NaN;

		for(int i=1; i<ntab; i++){
			deltaU = deltaU /con;
			
            a[0][i] = computeA(d, u, deltaU);
			
            fac = con2;
			
			for(int j=1; j<=i; j++){
				a[j][i] = (a[j-1][i]*fac - a[j-1][i-1])/(fac-1);
				fac = con2*fac;
				errt = Math.max(Math.abs(a[j][i]-a[j-1][i]), Math.abs(a[j][i]-a[j-1][i-1]));
				
				if (errt <= err){
					err = errt;
					d2fdu2 = a[j][i];
				}
			}
			
			if (Math.abs(a[i][i]-a[i-1][i-1]) >= safe*err){
				break;
			}
		}
		
		return d2fdu2;
	}

	public double computeA(int[] d, double[] u, double deltaU){
		
		if(d.length > 2){
			throw new RuntimeException("In computeA method; Derivative method can only do first and" +
					"second order derivative!!!");
		}
		
		if(d.length == 1){
			
			u[d[0]] =   deltaU;
			double fPlus = f(u);
			u[d[0]] = - deltaU;
			double fMinus= f(u);
			u[d[0]] = 0.0;
			
			setToInitialPosition();
			return (fPlus - fMinus)/(2.0*deltaU);
			
		} else {
			/*
			 *  SECOND PARTIAL DERIVATIVE
			 *  
			 *  f is the energy function
			 */
		
			if(d[0]==d[1]){
				/*
				 * 
				 *  when u_i = u_j
				 *  
				 *      d2f            f(u_i + deltaU) - 2*f(u_i) + f(u_i - deltaU) 
				 *  -----------  = ----------------------------------------------------- 
				 *   du_i*du_i                            deltaU^2
				 *  
				 */
				
				double f_init = f(u);
				u[d[0]] = deltaU;
				double f_up1 = f(u);
				u[d[0]] = - deltaU;
				double f_um1 = f(u);
				u[d[0]] = u[d[1]] = 0;
				
				setToInitialPosition();
				return (f_up1 - 2*f_init + f_um1)/(deltaU*deltaU);  
			
			} else {
				
				/*
				 *  when u_i, u_j
				 *  
				 *                  [  f(u_i+deltaU, u_j+deltaU) - f(u_i+deltaU, u_j-deltaU) - 
				 *     d2f                                    f(u_i-deltaU, u_j+deltaU) + f(u_i-deltaU, u_j-deltaU) ]
				 *  -----------  =  ------------------------------------------------------------------------------------ 
				 *   du_i*du_j                            4 * deltaU^2
				 *  
				 */
				
				u[d[0]] = deltaU;
				u[d[1]] = deltaU;
				double f_ip1jp1 = f(u);
				
				u[d[1]] = - deltaU;
				double f_ip1jm1 = f(u);
				
				u[d[0]] = - deltaU;
				double f_im1jm1 = f(u);
				
				u[d[1]] = deltaU;
				double f_im1jp1 = f(u);
				
				u[d[0]] = u[d[1]] = 0;
				
				setToInitialPosition();
				return (f_ip1jp1 - f_ip1jm1 - f_im1jp1 + f_im1jm1)/(4*deltaU*deltaU);
			}
		}
	}
	
	public void setToInitialPosition(){
		double newU[] = new double[coordinateDefinition.getCoordinateDim()];
		for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, newU);
		}
	}
	
	
	
	public double getFixedDeltaU() {
		return fixedDeltaU;
	}

	public void setFixedDeltaU(double fixedDeltaU) {
		this.fixedDeltaU = fixedDeltaU;
	}

	protected CoordinateDefinition coordinateDefinition;
	protected MeterPotentialEnergy meterPotential;
	protected Box box;
	protected double deltaU, errt, fac;
	protected double [][] a;
	double fixedDeltaU = 0.01;
	final int ntab = 10;
	final double con = 1.4;
	final double con2 = con*con;
	final double big = Double.MAX_VALUE;
	final double safe = 2.0;
	
	
}
