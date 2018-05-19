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
 * Finite Difference Derivative Class that is used by Conjugate Gradients related classes
 * 
 * This class returns the potential energy and the first derivative of lattice energy
 * PAY ATTENTION to variable: fixedDeltaU 
 * 
 * @author taitan
 *
 */
public class FiniteDifferenceDerivativeCGNitrogenBeta{
	
	public FiniteDifferenceDerivativeCGNitrogenBeta(Box box, PotentialMaster potentialMaster,
                                                    CoordinateDefinition coordinateDefinition){
		this.coordinateDefinition = coordinateDefinition;
		this.box = box;
		meterPotential = new MeterPotentialEnergy(potentialMaster, box);
		
		a = new double[ntab][ntab];
	}
 	
	public double f(double[] parameter) {
		double[] newU = generateNewU(parameter); 
		for (int cell=0; cell<coordinateDefinition.getBasisCells().length; cell++){
			IMoleculeList molecules = coordinateDefinition.getBasisCells()[cell].molecules;
			coordinateDefinition.setToU(molecules, newU);
		}
		
		return meterPotential.getDataAsScalar();
	}

	public double dphi_du(int[] d, double parameter[]) {
		
		/*
		 * d[i] = the i-th elements of the derivative
		 * u[i] = the generalized coordinate/ relative displacement from 
		 * 			nominal position 
		 */
        deltaU = fixedDeltaU;

        a[0][0] = computeA(d, parameter, deltaU);
     
		double err = big;
        double dfdu = Double.NaN;
		
		for(int i=1; i<ntab; i++){
			deltaU = deltaU /con;
            a[0][i] = computeA(d, parameter, deltaU);
            
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

	public double computeA(int[] d, double[] parameter, double deltaU){
		
//		if(d.length > 2){
//			throw new RuntimeException("In computeA method; Derivative method can only do first and" +
//					"second order derivative!!!");
//		}
		
		int k=Integer.MAX_VALUE;
		for (int i=0; i<d.length; i++){
			if(d[i]==1){
				k=i;
				break;
			}
		}
		double initParami = parameter[k];
		
		parameter[k] +=  deltaU;
		double fPlus = f(parameter);
		parameter[k] -=  deltaU;
		double fMinus= f(parameter);
		parameter[k] = initParami;
		
		setToInitialPosition(parameter);
		return (fPlus - fMinus)/(2.0*deltaU);
			
	}
	
	public void setToInitialPosition(double[] parameter){
		double[] newU = generateNewU(parameter);
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

	public double[] generateNewU (double[] parameter){
		int numDOF = coordinateDefinition.getCoordinateDim();
		int dofPerMol = numDOF/box.getMoleculeList().size();
		
		int nC = (int)Math.pow(box.getMoleculeList().size()/1.9999999999, 1.0/3.0);
		int nCelldofinZ = nC*2*dofPerMol;
		double[] newU = new double[numDOF];
		
		boolean isUnitCellA = true;
		int counter = 1;
		for(int i=0; i<newU.length; i+=10){
			
			if(i>0 && i%nCelldofinZ == 0){
				isUnitCellA = !isUnitCellA;
				
				if((nC%2 == 1) && (i==nCelldofinZ*nC*counter)){
					isUnitCellA = !isUnitCellA;
					++ counter;
				}
			}
			
			if(isUnitCellA){
				newU[i] = parameter[0];
				newU[i+1] = parameter[1];
				newU[i+2] = parameter[2];
				newU[i+3] = parameter[3];
				newU[i+4] = parameter[4];
			
				newU[i+5] = parameter[5];
				newU[i+6] = parameter[6];
				newU[i+7] = parameter[7];
				newU[i+8] = parameter[8];
				newU[i+9] = parameter[9];
			
			} else {
				newU[i] = parameter[10];
				newU[i+1] = parameter[11];
				newU[i+2] = parameter[12];
				newU[i+3] = parameter[13];
				newU[i+4] = parameter[14];
			
				newU[i+5] = parameter[15];
				newU[i+6] = parameter[16];
				newU[i+7] = parameter[17];
				newU[i+8] = parameter[18];
				newU[i+9] = parameter[19];
			}
		}
		return newU;
	}
	
	
	
	protected CoordinateDefinition coordinateDefinition;
	protected MeterPotentialEnergy meterPotential;
	protected Box box;
	protected double deltaU, errt, fac;
	protected double [][] a;
	double fixedDeltaU = 1e-10;
	final int ntab = 10;
	final double con = 1.4;
	final double con2 = con*con;
	final double big = Double.MAX_VALUE;
	final double safe = 2.0;
	
	
}
