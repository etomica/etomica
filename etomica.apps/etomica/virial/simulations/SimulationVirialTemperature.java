package etomica.virial.simulations;

import etomica.virial.*;
import etomica.*;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class SimulationVirialTemperature extends SimulationVirial {

	/**
	 * Constructor for SimulationVirialTemperature.
	 * @param nMolecules
	 * @param simTemperature
	 */
	public SimulationVirialTemperature(int nMolecules, double simTemperature) {
		super(nMolecules, simTemperature);
	}

	/**
	 * Constructor for SimulationVirialTemperature.
	 * @param space
	 * @param nMolecules
	 * @param simTemperature
	 */
	public SimulationVirialTemperature(
		Space space,
		int nMolecules,
		double simTemperature) {
		super(space, nMolecules, simTemperature);
	}

	public static void main(String[] args) {
		double targetTemperature = 1.3;  //temperature of calculated virial coefficient
		double simTemperature = 0.9*targetTemperature; //temperature governing sampling of configurations
		double sigmaHSMod = sigmaLJ1B(1.0/simTemperature); //range in which modified-f for sampling will apply abs() function
		int nMolecules = 3;
		SimulationVirialTemperature sim = new SimulationVirialTemperature(nMolecules, simTemperature);		
		sim.species().setDiameter(sigmaHSMod);
		sim.integrator().setDoSleep(true);
		
		//set up simulation potential
		P0Cluster p2 = new P0Cluster(sim.hamiltonian.potential, sim.pairs());
		p2.setTemperature(simTemperature);		
		P2LennardJones p2LJ = new P2LennardJones(p2);
		MayerModified fSim = new MayerModified(p2LJ, sigmaHSMod);
		Cluster ringSim = new etomica.virial.cluster.Ring(nMolecules, 1.0, fSim);
		p2.setCluster(ringSim);				
		sim.setSimPotential(p2);
		
		MayerGeneral f = new MayerGeneral(p2LJ);
		Cluster ring = new etomica.virial.cluster.Ring(nMolecules, 1.0, f);
		Cluster[] clusters = new Cluster[] {ring};
		Cluster refCluster = ring;
		double refIntegral = 1.0;
		double refTemperature = simTemperature;
		
		MeterVirial meterVirial = new MeterVirial(sim, sim.pairs(), refTemperature, refCluster, refIntegral, 
													clusters, p2);
		meterVirial.setTemperature(targetTemperature);
		sim.setMeterVirial(meterVirial);

	}
}
