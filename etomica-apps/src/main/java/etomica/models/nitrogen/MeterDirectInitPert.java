/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.util.random.IRandom;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.integrator.IntegratorBox;
import etomica.units.Null;

/**
 * Meter used for overlap sampling in the target-sampled system.  The meter
 * measures the ratio of the Boltzmann factors for the reference and target
 * potentials.
 * 
 * 
 * 
 * @author Tai Boon Tan
 */
public class MeterDirectInitPert extends DataSourceScalar {
	
	public MeterDirectInitPert(IntegratorBox integrator, MeterPotentialEnergy meterPotentialEnergy, CoordinateDefinitionNitrogen coordinateDef, IRandom random) {
    	super("Scaled Unit",Null.DIMENSION);
         
    	meterEnergy = new MeterPotentialEnergyFromIntegrator(integrator);
        this.integrator = integrator;
        this.meterPotentialEnergy = meterPotentialEnergy;
        this.coordinateDef = coordinateDef;
        this.random = random;
        
		newU = new double[coordinateDef.getCoordinateDim()];
		initU = new double[coordinateDef.getCoordinateDim()];

    }
    
	public double getDataAsScalar() {
		double uTransOnly = meterEnergy.getDataAsScalar();
		
		initU = coordinateDef.calcU(meterPotentialEnergy.getBox().getMoleculeList());
		
		/*
		 * Generating random u3 and u4 values that satisfy the u_max value (constraint angle)
		 * 
		 * 	u3^2 + u4^2 = 2[ 1- cos(theta) ]
		 *  at small theta limit, the equation becomes:
		 *  u3^2 + u4^2 = theta^2
		 *  
		 *  u3^2 + u4^2 = u_max^2
		 *  u_max^2 = 2[ 1 - cos(theta) ] = theta^2
		 *  u_max = sqrt(2[1-cos(theta)]) = theta
		 *  
		 */
		for(int i=0; i<coordinateDef.getCoordinateDim(); i+=5){
			for(int j=0; j<3; j++){
				newU[i+j] = initU[i+j];
			}
			
			double u3 = 2.0;
			double u4 = 2.0;
			while ( (u3*u3+u4*u4) > 2*(1-Math.cos(constraintAngle))){
				u3 = (2*random.nextDouble() - 1.0)*(constraintAngle/180.0) *Math.PI;
				u4 = (2*random.nextDouble() - 1.0)*(constraintAngle/180.0) *Math.PI;
			} 
			newU[i+3] = u3;
			newU[i+4] = u4;
			
		}
		coordinateDef.setToU(meterPotentialEnergy.getBox().getMoleculeList(), newU);
		
		double uTransAndRotate = meterPotentialEnergy.getDataAsScalar();
    	coordinateDef.setToU(meterPotentialEnergy.getBox().getMoleculeList(), initU);
		
    	if(Double.isNaN(uTransOnly) || Double.isNaN(uTransAndRotate)){
    		throw new RuntimeException("<MeterDirectInitPert> energy is NaN!!!!!!!!!!!!");
    	}
    	 
    	//System.out.println(uSampled + " " + uMeasured + " "+ (uMeasured-uSampled) + " " + Math.exp(-(uMeasured - uSampled) / integrator.getTemperature()));
    	return Math.exp(-(uTransAndRotate - uTransOnly) / integrator.getTemperature());
	}
	

	public double getConstraintAngle() {
		return constraintAngle;
	}

	public void setConstraintAngle(double constraintAngle) {
		this.constraintAngle = constraintAngle;
	}
	
	private static final long serialVersionUID = 1L;
    protected final MeterPotentialEnergyFromIntegrator meterEnergy;
    protected final MeterPotentialEnergy meterPotentialEnergy;
    protected final IntegratorBox integrator;
    protected double constraintAngle;
    protected CoordinateDefinitionNitrogen coordinateDef;
    protected IRandom random;
    protected double[] newU, initU;

}
