package etomica.virial.overlap;

import etomica.*;
import etomica.space3d.Boundary;
import etomica.virial.*;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;

/**
 * @author kofke
 *
 * Overlap-sampling simulation of the target system.
 */
public class SimulationOSTarget extends SimulationGraphic {

	/**
	 * Constructor for SimulationOverlap.
	 */
	public SimulationOSTarget(Space space, double temperature, boolean targetPositive,
			Cluster refCluster, Cluster targetCluster) {
		super(space);

		int nMolecules = targetCluster.pointCount();
		simTemperature = temperature;
		double sigmaHSRef = 1.0*((MayerHardSphere)refCluster.bondGroup()[0].f).getSigma();
		boolean refPositive = !refCluster.hasOddBondCount();

		phase = new PhaseCluster(this);
		phase.setBoundary(space.makeBoundary(Boundary.NONE));	
		species = new SpeciesSpheresMono(this);
		species.setNMolecules(nMolecules);
		species.setDiameter(sigmaHSRef);
		this.elementCoordinator.go();
		
		Controller controller = new Controller(this);		
		DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(this);
		integrator = new IntegratorMC(this);
		integrator.setSleepPeriod(1);
		integrator.setTemperature(simTemperature);
		MCMoveAtom mcMoveAtom1 = new MeterVirial.MyMCMoveAtom(integrator);
		for(int n=2; n<nMolecules; n++) {
			new MCMoveAtomMulti(integrator, n);
		}
		
		//set up simulation potential for reference cluster
		P0ClusterSigned p0 = new P0ClusterSigned(this.hamiltonian.potential, targetCluster);
		p0.setSignPositive(targetPositive);
	
	    ConfigurationCluster configuration = new ConfigurationCluster(this);
	    configuration.setPhase(phase);
	    configuration.setCluster(targetCluster);
	    configuration.setSignPositive(targetPositive);
	    phase.setConfiguration(configuration);						
			
	    meter = new MeterOverlapTarget(this, refPositive, 
	    										targetCluster, refCluster);
	    meter.setTemperature(simTemperature);
		meter.setLabel("Target "+(targetPositive?"Positive":"Negative"));
		
	    DisplayPlot clusterPlot = new DisplayPlot(this);
	    clusterPlot.setDataSources(meter);
	    clusterPlot.setWhichValue(MeterAbstract.AVERAGE);
		
	    this.elementCoordinator.go();
	    clusterPlot.setDataSources(meter);
	    clusterPlot.setLabel("Target "+ (targetPositive?"Positive":"Negative"));

		DisplayPhase display = new DisplayPhase();
		ColorSchemeByType.setColor(species, java.awt.Color.green);
		
		this.elementCoordinator.go();
	}//end of constructor
	
	private MeterOverlapTarget meter;
	private double simTemperature;
	private double refTemperature;
	private P0Cluster p2;
	private SpeciesSpheresMono species;
	protected IntegratorMC integrator;
	private PhaseCluster phase;

	public Phase phase() {return phase;}
	
	public IntegratorMC integrator() {return integrator;}
	/**
	 * Returns the meterVirial.
	 * @return MeterVirial
	 */
	public MeterOverlapTarget getMeter() {
		return meter;
	}

	/**
	 * Returns the refTemperature.
	 * @return double
	 */
	public double getRefTemperature() {
		return refTemperature;
	}

	/**
	 * Returns the simTemperature.
	 * @return double
	 */
	public double getSimTemperature() {
		return simTemperature;
	}

	/**
	 * Sets the refTemperature.
	 * @param refTemperature The refTemperature to set
	 */
	public void setRefTemperature(double refTemperature) {
		this.refTemperature = refTemperature;
	}

	/**
	 * Sets the simTemperature.
	 * @param simTemperature The simTemperature to set
	 */
	public void setSimTemperature(double simTemperature) {
		this.simTemperature = simTemperature;
	}
	
	public P0Cluster getSimPotential() {
		return p2;
	}
	
	public void setSimPotential(P0Cluster p2) {
		this.p2 = p2;
	}
	
	public SpeciesSpheresMono species() {
		return species;
	}

}
