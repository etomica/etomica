package etomica.virial.overlap;

import etomica.*;
import etomica.graphics.*;
import etomica.virial.*;

/**
 * @author kofke
 *
 * Simulation implementing the overlap-sampling approach to evaluating a cluster
 * diagram.
 */
public class SimulationOverlap extends SimulationGraphic {

	
	public SimulationOverlap(Space space, double temperature, 
								Cluster targetCluster, Cluster refCluster) {
		super(space);

		Default.TRUNCATE_POTENTIALS = false;		
		
		int nMolecules = targetCluster.pointCount();
		simTemperature = temperature;
		double sigmaHSRef = 1.0*((MayerHardSphere)refCluster.bondGroup()[0].f).getSigma();
		boolean refPositive = !refCluster.hasOddBondCount();

///////// reference-system simulation				
		this = this;
		phase = new Phase(this);
		phase.setBoundary(space.makeBoundary(Space3D.Boundary.NONE));	
		species = new SpeciesSpheresMono(this);
		species.setNMolecules(nMolecules);
		species.setDiameter(sigmaHSRef);
		this.elementCoordinator.go();
		pairs = new PairSet(((AtomTreeNodeGroup)phase.getAgent(species).node).childList);

		refCluster.setPairSet(pairs);
		targetCluster.setPairSet(pairs);		
		
		Controller controller = new Controller(this);		
		DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(this);
		integrator = new IntegratorMC(this);
		integrator.setSleepPeriod(1);
		integrator.setTemperature(simTemperature);
		MCMoveAtom mcMoveAtom1 = new MeterVirial.MyMCMoveAtom(integrator);
		MCMoveAtomMulti mcMoveAtom2 = new MCMoveAtomMulti(integrator,2);
		for(int n=3; n<nMolecules; n++) {
			new MCMoveAtomMulti(integrator, n);
		}
		
		//set up simulation potential for reference cluster
		P2ClusterSigned p2 = new P2ClusterSigned(this.hamiltonian.potential, pairs);
		p2.setCluster(refCluster);
		p2.setSignPositive(refPositive);
		p2.setTemperature(simTemperature);			

	  boolean simulatingTarget = false;
	  boolean targetPositive = false;
	  Cluster simCluster = simulatingTarget ? targetCluster : refCluster;
	  Cluster nonSimCluster = simulatingTarget ? refCluster : targetCluster;
	
	  ConfigurationCluster configuration = new ConfigurationCluster(this);
	  configuration.setPhase(phase);
	  configuration.setCluster(simCluster);
	  configuration.setSignPositive(simulatingTarget ? targetPositive : refPositive);
	  phase.setConfiguration(configuration);						
			
	  MeterOverlapReference meter = new MeterOverlapReference(this, simCluster, nonSimCluster);
	  meter.setTemperature(simTemperature);
	  meter.setActive(true);
		
	  DisplayPlot clusterPlot = new DisplayPlot(this);
	  clusterPlot.setDataSources(meter.allMeters());
	  clusterPlot.setWhichValue(MeterAbstract.AVERAGE);
		
	  this.elementCoordinator.go();
	  clusterPlot.setDataSources(meter.allMeters());
	  clusterPlot.setLabel("Reference");

		DisplayPhase display = new DisplayPhase();
		ColorSchemeByType.setColor(species, java.awt.Color.green);
		
		this.elementCoordinator.go();
		
	}
	
	private MeterOverlap meter;
	private double simTemperature;
	private double refTemperature;
	private PairSet pairs;
	private P2Cluster p2;
	private SpeciesSpheresMono species;
	protected IntegratorMC integrator;
	private Phase phase;
	
	public Phase phase() {return phase;}
	
	public IntegratorMC integrator() {return integrator;}
	/**
	 * Returns the meterVirial.
	 * @return MeterVirial
	 */
	public MeterOverlap getMeter() {
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
	
	public PairSet pairs() {
		return pairs;
	}
		
	public P2Cluster getSimPotential() {
		return p2;
	}
	
	public void setSimPotential(P2Cluster p2) {
		this.p2 = p2;
	}
	
	public SpeciesSpheresMono species() {
		return species;
	}
		
}
