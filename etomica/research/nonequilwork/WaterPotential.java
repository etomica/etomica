package etomica.research.nonequilwork;

import etomica.models.water.*;
import etomica.*;
import etomica.units.Kelvin;
import etomica.graphics.*;
//import etomica.log.*;

public class WaterPotential extends SimulationGraphic {
    public WaterPotential() {
        super(new Space3D());
        Simulation.instance = this;
        Default.ATOM_SIZE = 3.1670;
        DefaultGraphic.ATOM_COLOR = java.awt.Color.RED;
        int nMolecules = 32;
        double density = 6.02/180 ;
        PotentialTruncationSimple pts = new PotentialTruncationSimple(this.space, 6.124);
	    
	    IntegratorMC integrator = new IntegratorMC(this);
	    integrator.setTemperature(Kelvin.UNIT.toSim(298.));
	    MCMoveMolecule mcMoveMolecule = new MCMoveMolecule(integrator);
	    MCMoveRotateMolecule3D mcRotateMolecule = new MCMoveRotateMolecule3D(integrator);
	    SpeciesWater water = new SpeciesWater(this);
	    water.setNMolecules(nMolecules);
	    
	    PotentialOSInsert potentialOS = new PotentialOSInsert(this.hamiltonian.potential);
		P2WaterSPC p2Water = new P2WaterSPC(potentialOS, pts);
		potentialOS.setPotential(p2Water);
		potentialOS.setSpecies(water);
	    
	    Phase phase = new Phase(this); 
        phase.setLrcEnabled(false);  
		p2Water.setBoundary((Space3D.Boundary)phase.getBoundary());
        elementCoordinator.go();
        
		potentialOS.setTestMolecule(phase.molecule(0));
		      
        phase.getCentralImageEnforcer().setApplyToMolecules(true);
	    phase.setDensity(density*0.97);
	    integrator.setInterval(nMolecules);
	    integrator.setDoSleep(true);
	    
	    Controller controller = new Controller(this);
        MeterCycles meterCycles = new MeterCycles();
	    
		MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(this);
		meterEnergy.setActive(true);
	
//	    TableDatumSources table = new TableDatumSources();
//        table.setDatumSources(new MeterScalar[]{meterCycles, meterEnergy});
//        table.setWhichValues(MeterAbstract.AVERAGE);
//        LoggerAbstract logger = new LoggerAbstract(table);
//        logger.setFileName(name);
//        logger.setFileNameSuffix(".xls");
//        logger.setUpdateInterval(50);
//        logger.setPrecision(10);
//        logger.setCloseFileEachTime(true);
//

		DisplayPhase display = new DisplayPhase();
		display.setPhase(phase);
		
		elementCoordinator.go();
		
		this.makeAndDisplayFrame();
		
//		controller.start();
//		controller.setMaxSteps(900000000);
    }
    
  
    public static void main(String[] args) {
        WaterPotential sim = new WaterPotential();
    }//end of main
    
}
