package etomica.virial.simulations;

import etomica.*;
import etomica.graphics.*;
import etomica.virial.*;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class SimulationVirial extends SimulationGraphic {

	/**
	 * Default constructor, using a 3D space.
	 * @see java.lang.Object#Object()
	 */
	public SimulationVirial(int nMolecules, double simTemperature) {
		this(new Space3D(), nMolecules, simTemperature);
	}
	
	/**
	 * Constructor for SimulationVirial.
	 */
	public SimulationVirial(Space space, int nMolecules, double simTemperature) {
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
	private MeterVirial meterVirial;
	private double simTemperature;
	private double refTemperature;
	private PairSet pairs;
	private P2Cluster p2;
	private SpeciesSpheresMono species;
	protected IntegratorMC integrator;
	
	public IntegratorMC integrator() {return integrator;}
	/**
	 * Returns the meterVirial.
	 * @return MeterVirial
	 */
	public MeterVirial getMeterVirial() {
		return meterVirial;
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
	public void setMeterVirial(MeterVirial meterVirial) {
		this.meterVirial = meterVirial;
		meterVirial.setPairs(pairs);

		MeterDatumSourceWrapper bMeter = new MeterDatumSourceWrapper(meterVirial);
		bMeter.setHistorying(true);
		DisplayPlot bPlot = new DisplayPlot(this);
		bPlot.setDataSources(bMeter.getHistory());
		bPlot.setWhichValue(MeterAbstract.CURRENT);
		bMeter.getHistory().setNValues(1000);
		bPlot.setLabel("B running average");
		
		DisplayPlot clusterPlot = new DisplayPlot(this);
		meterVirial.setHistorying(true);
		MeterDatumSourceWrapper[] clusterMeter = new MeterDatumSourceWrapper[4];
		for(int i=0; i<2; i++) {
			clusterMeter[i] = new MeterDatumSourceWrapper(meterVirial.allMeters()[i]);
			clusterMeter[i].setLabel("Meter"+i);
			clusterMeter[i].setHistorying(true);
			clusterMeter[i].getHistory().setNValues(1000);
		}
		clusterPlot.setDataSources(new DataSource[] {
									clusterMeter[0].getHistory(),
									clusterMeter[1].getHistory()
									});
//		clusterPlot.setDataSources(meterVirial.getDataSources("History"));
//		clusterPlot.setWhichValue(MeterAbstract.AVERAGE);
		clusterPlot.setLabel("Cluster integrals");
		
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
		
	public P2Cluster getSimPotential() {
		return p2;
	}
	
	public void setSimPotential(P2Cluster p2) {
		this.p2 = p2;
	}
	
	public SpeciesSpheresMono species() {
		return species;
	}
		
	public static void main(String[] args) {
		double temperature = 1.3;  //temperature of calculated virial coefficient
		double simTemperature = 1.0*temperature; //temperature governing sampling of configurations
		double sigmaHSRef = sigmaLJ1B(1.0/simTemperature);  //diameter of reference HS system
		double sigmaHSMod = sigmaLJ1B(1.0/simTemperature); //range in which modified-f for sampling will apply abs() function
		System.out.println((1./simTemperature)+" "+sigmaHSRef);
		System.out.println((1./temperature)+" "+sigmaHSMod);
		int nMolecules = 3;
//		int nMolecules = 4;
		SimulationVirial sim = new SimulationVirial(nMolecules, simTemperature);
		
		sim.species().setDiameter(sigmaHSRef);
		
		//set up simulation potential
		P2Cluster p2 = new P2Cluster(sim.hamiltonian.potential, sim.pairs());
		p2.setTemperature(simTemperature);		
		P2LennardJones p2LJ = new P2LennardJones(p2);
		MayerModified f = new MayerModified(p2LJ, sigmaHSMod);
		Cluster ring = new etomica.virial.cluster.Ring(nMolecules, 1.0, f);
		p2.setCluster(ring);				
		sim.setSimPotential(p2);
		
		MeterVirial meterVirial = (nMolecules==3) ? (MeterVirial)new MeterVirialB3(sim, sim.pairs(), sigmaHSRef, p2, p2LJ)
												  : (MeterVirial)new MeterVirialB4(sim, sim.pairs(), sigmaHSRef, p2, p2LJ);
		meterVirial.setTemperature(temperature);
		sim.setMeterVirial(meterVirial);
	}//end of main
}
