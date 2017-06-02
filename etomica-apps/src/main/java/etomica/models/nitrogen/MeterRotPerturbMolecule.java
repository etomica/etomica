/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.api.ISpecies;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.integrator.IntegratorMC;
import etomica.space.Space;
import etomica.units.Null;

/**
 * Meter that returns a Boltzmann factor difference between 2 systems: system 1 and 0
 * 1: (Target System) sampled with an angle constraint, e.g. one degree angle
 * 0: (Reference System) with no rotational d.o.f.
 *  
 * @author Tai Boon Tan
 *
 */
public class MeterRotPerturbMolecule extends DataSourceScalar {

	private static final long serialVersionUID = 1L;
	protected final MeterPotentialEnergy meterPotentialMeasured;
	protected final MeterPotentialEnergyFromIntegrator meterPotentialSampled;
    protected final Box secondaryBox;
    protected double latticeEnergy;
    protected CoordinateDefinitionNitrogen primaryCoordDef, secondaryCoordDef;
    
    public MeterRotPerturbMolecule(IntegratorMC integrator, PotentialMaster potentialMaster, ISpecies species, Space space, Simulation sim, CoordinateDefinitionNitrogen coordinateDef) {
        super("Scaled Energy unit", Null.DIMENSION);
        this.primaryCoordDef = coordinateDef;
        
        Box realBox = coordinateDef.getBox();
        secondaryBox = new Box(space);
        sim.addBox(secondaryBox);
       
        secondaryBox.setNMolecules(species, realBox.getNMolecules(species));
        secondaryBox.setBoundary(realBox.getBoundary());
     
        secondaryCoordDef = new CoordinateDefinitionNitrogen(sim, secondaryBox, coordinateDef.getPrimitive(), coordinateDef.getBasis(), space);
        secondaryCoordDef.setIsBeta();
        secondaryCoordDef.setOrientationVectorBeta(space);
        secondaryCoordDef.initializeCoordinates(new int[]{1,1,1});
        
        meterPotentialMeasured = new MeterPotentialEnergy(potentialMaster);
        meterPotentialSampled = new MeterPotentialEnergyFromIntegrator(integrator);
        
    }

	public double getDataAsScalar() {
		
		double[] sampledCoord = primaryCoordDef.calcU(primaryCoordDef.getBox().getMoleculeList());
		double[] transCoord = new double[sampledCoord.length];
				
		for(int i=0; i<sampledCoord.length; i++){
			if(i>0 && (i%5==3 || i%5==4)){
				transCoord[i] = 0.0;
				
			} else{
				transCoord[i] = sampledCoord[i];
				
			}
		}
		
		secondaryCoordDef.setToU(secondaryBox.getMoleculeList(), transCoord);
		meterPotentialMeasured.setBox(secondaryBox);
		
		double sampledEnergy = meterPotentialSampled.getDataAsScalar();
		double measuredEnergy = meterPotentialMeasured.getDataAsScalar();
		
		double chi = Math.exp(-(measuredEnergy-sampledEnergy)/meterPotentialSampled.getIntegrator().getTemperature()); 
		//System.out.println(sampledEnergy+" "+measuredEnergy + " "+(measuredEnergy-sampledEnergy) + " " + chi);
		return chi;
	}
	
	public void setLatticeEnergy(double latticeEnergy){
		this.latticeEnergy = latticeEnergy;
	}

}
