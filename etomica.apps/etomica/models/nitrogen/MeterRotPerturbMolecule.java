package etomica.models.nitrogen;

import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.integrator.IntegratorMC;
import etomica.space.ISpace;
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
    protected double temperature;
    protected double[] otherTemperatures;
    protected final IBox secondaryBox;
    protected CoordinateDefinitionNitrogen primaryCoordDef, secondaryCoordDef;
    
    public MeterRotPerturbMolecule(IntegratorMC integrator, IPotentialMaster potentialMaster, ISpecies species, ISpace space, ISimulation sim, CoordinateDefinitionNitrogen coordinateDef) {
        super("Scaled Energy unit", Null.DIMENSION);
        this.primaryCoordDef = coordinateDef;
        
        IBox realBox = coordinateDef.getBox();
        secondaryBox = new Box(space);
        secondaryBox.setBoundary(realBox.getBoundary());
        secondaryBox.setNMolecules(species, realBox.getMoleculeList().getMoleculeCount());
        sim.addBox(secondaryBox);
        
        secondaryCoordDef = new CoordinateDefinitionNitrogen(sim, secondaryBox, coordinateDef.getPrimitive(), coordinateDef.getBasis(), space);
        secondaryCoordDef.setIsBeta();
        secondaryCoordDef.setOrientationVectorBeta(space);
        secondaryCoordDef.initializeCoordinates(new int[]{1,1,1});
        
        meterPotentialMeasured = new MeterPotentialEnergy(potentialMaster);
        meterPotentialMeasured.setBox(secondaryBox);
        
        meterPotentialSampled = new MeterPotentialEnergyFromIntegrator(integrator);
        
    }

	public double getDataAsScalar() {
		
		double[] sampledCoord = primaryCoordDef.calcU(primaryCoordDef.getBox().getMoleculeList());
		double[] transCoord = new double[sampledCoord.length];
		       
		for(int i=0; i<sampledCoord.length; i++){
			if(i==0 && (i%3==0 || i%4==0)){
				transCoord[i] = 0.0;
				
			}else{
				transCoord[i] = sampledCoord[i];
				
			}
		}
		
		secondaryCoordDef.setToU(secondaryBox.getMoleculeList(), transCoord);
		
		double sampledEnergy = meterPotentialSampled.getDataAsScalar();
		double measuredEnergy = meterPotentialMeasured.getDataAsScalar();
		
		double chi = Math.exp(-(measuredEnergy-sampledEnergy)/meterPotentialSampled.getIntegrator().getTemperature()); 
		
		return chi;
	}

}
