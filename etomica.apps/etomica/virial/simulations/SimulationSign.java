package etomica.virial.simulations;

import etomica.*;
import etomica.graphics.*;
import etomica.virial.*;
import etomica.virial.cluster.*;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class SimulationSign extends SimulationGraphic {

	/**
	 * Default constructor, using a 3D space.
	 * @see java.lang.Object#Object()
	 */
	public SimulationSign(int nMolecules, double simTemperature) {
		this(new Space3D(), nMolecules, simTemperature);
	}
	
	/**
	 * Constructor for SimulationVirial.
	 */
	public SimulationSign(Space space, int nMolecules, double simTemperature) {
		super(space);

		Default.makeLJDefaults();
		Default.TRUNCATE_POTENTIALS = false;
		
		final Phase phase = new Phase(this);
		phase.setBoundary(this.space.makeBoundary(Space3D.Boundary.NONE));	
		species = new SpeciesSpheresMono(this);
		species.setNMolecules(nMolecules);
		elementCoordinator.go();
		pairs = new PairSet(((AtomTreeNodeGroup)phase.getAgent(species).node).childList);
		
		Controller controller = new Controller(this);		
		DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(this);
		integrator = new IntegratorMC(this);
		integrator.setSleepPeriod(1);
		integrator.setDoSleep(false);
		integrator.setTemperature(simTemperature);
		MCMoveAtom mcMoveAtom = new MeterVirial.MyMCMoveAtom(integrator);
		
		DisplayPhase display = new DisplayPhase();
		ColorSchemeByType.setColor(species, java.awt.Color.green);
		
		this.elementCoordinator.go();
		
		Space3D.Vector origin = new Space3D.Vector(5.,5.,5.);
		SpeciesAgent speciesAgent = phase.getAgent(species);
		speciesAgent.coord.translateTo(new Space3D.Vector(5.,5.,5.));
		speciesAgent.firstMolecule().coord.translateTo(new Space3D.Vector(5.,5.,5.));
				
		AtomList childList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
		AtomIteratorList list1 = new AtomIteratorList(childList);
		list1.reset();
		while(list1.hasNext()) list1.next().coord.translateTo(origin);
	}
	
	public static double sigmaLJ1 = 1.28412293285;//  ( 4 + 4 sqrt(1-Ln(2)) ) / Ln(4))^(1/6), which is where f(r) = 1 for LJ

	/**
	 * Returns the separation r for the LJ potential at which beta*u(r) = -Ln(2)
	 * (so that f(r) = 1)
	 * @param beta
	 * @return double
	 */
	public static double sigmaLJ1B(double beta) {
		double log2 = Math.log(2.0);
		if(beta <= log2) return Math.pow(2.0, 1./6.);
		else return Math.pow( 2.0*(beta + Math.sqrt(beta*(beta-log2)) )/log2, 1./6.);
	}
	private MeterSign meterSign;
	private double simTemperature;
	private double refTemperature;
	private PairSet pairs;
	private P0Cluster p2;
	private SpeciesSpheresMono species;
	protected IntegratorMC integrator;
	
	public IntegratorMC integrator() {return integrator;}
	/**
	 * Returns the meterVirial.
	 * @return MeterVirial
	 */
	public MeterSign getMeterSign() {
		return meterSign;
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
	 * Sets the meterVirial and configures displays for it.
	 * @param meterVirial The meterVirial to set
	 */
	public void setMeterSign(MeterSign meterSign) {
		this.meterSign = meterSign;
		meterSign.setPairSet(pairs);

		MeterDatumSourceWrapper bMeter = new MeterDatumSourceWrapper(meterSign);
		bMeter.setWhichValue(MeterAbstract.AVERAGE);
		bMeter.setHistorying(true);
		meterSign.setHistorying(true);
		DisplayPlot bPlot = new DisplayPlot(this);
//		bPlot.setDataSources(new DataSource[] {bMeter.getHistory(), meterSign.getHistory()});
		bPlot.setDataSources(new DataSource[] {bMeter.getHistory()});
//		bPlot.setWhichValue(MeterAbstract.AVERAGE);
		meterSign.getHistory().setNValues(1000);
		bMeter.getHistory().setNValues(1000);
		bPlot.setLabel("Sign running average");
				
		this.elementCoordinator.go();
		this.makeAndDisplayFrame();
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
		
	public P0Cluster getSimPotential() {
		return p2;
	}
	
	public void setSimPotential(P0Cluster p2) {
		this.p2 = p2;
	}
	
	public SpeciesSpheresMono species() {
		return species;
	}
		
	public static void main(String[] args) {
		double simTemperature = 1.3; //temperature governing sampling of configurations
		double sigmaHSMod = 0.0;//sigmaLJ1B(1.0/simTemperature); //range in which modified-f for sampling will apply abs() function
		int nMolecules = 4;
//		int nMolecules = 4;
		SimulationSign sim = new SimulationSign(nMolecules, simTemperature);
		
		sim.species().setDiameter(1.0);
//		sim.integrator().setDoSleep(true);
		
		//set up simulation potential
		P0Cluster p2 = new P0Cluster(sim.hamiltonian.potential, sim.pairs());
		p2.setTemperature(simTemperature);		
		P2LennardJones p2LJ = new P2LennardJones(p2);
		MayerModified f = new MayerModified(p2LJ, sigmaHSMod);
//		Cluster ring = new etomica.virial.cluster.Full(nMolecules, 1.0, f);
		Cluster ring = new ReeHoover(nMolecules, 1.0, new Cluster.BondGroup(f,Standard.D4));
		p2.setCluster(ring);				
		sim.setSimPotential(p2);
		
		MayerFunction fLJ = new MayerGeneral(p2LJ);
//		Cluster cluster = new etomica.virial.cluster.Full(nMolecules, 1.0, fLJ);
		Cluster cluster = new ReeHoover(nMolecules, 1.0, new Cluster.BondGroup(fLJ,Standard.D4));
		MeterSign meterSign = new MeterSign(sim, sim.pairs(), cluster);
		sim.setMeterSign(meterSign);
	}//end of main
}
