package etomica.virial;

import etomica.*;
import etomica.virial.cluster.Ring;
import etomica.graphics.*;

/**
 * @author David Kofke
 *
 */
public class MeterVirialB3 extends MeterVirial {

	/**
	 * Constructor for MeterVirialB3.
	 * @param sim
	 * @param pairSet
	 * @param sigma
	 * @param simulationPotential
	 * @param targetPotential
	 */
	public MeterVirialB3(Simulation sim, PairSet pairSet, double refSigma, P2Cluster simulationPotential, Potential2 targetPotential) {
		super(
			sim,
			pairSet,
			B3HS(refSigma),
			new Ring(3, -1./3., new MayerHardSphere(refSigma)),
			B3Clusters(targetPotential),
			simulationPotential);
	}
	
	public static double B3HS(double sigma) {
		double b0 = 2.*Math.PI/3 * sigma*sigma*sigma;
		return -5./8. * b0 * b0;
	}
	
	public static Cluster[] B3Clusters(Potential2 potential) {
		return new Cluster[] {new Ring(3, -1./3., new MayerGeneral(potential))};
	}

	public static void main(String[] args) {
	
		Default.makeLJDefaults();
		Default.TRUNCATE_POTENTIALS = false;
		int nMolecules = 3;
		double sigmaLJ1 = 1.28412293285;//  ( 4 + 4 sqrt(1-Ln(2)) ) / Ln(4))^(1/6), which is where f(r) = 1 for LJ
		double sigmaHSRef = sigmaLJ1;  //diameter of reference HS system
		double sigmaHSMod = sigmaLJ1; //range in which modified-f for sampling will apply abs() function
		double temperature = 1.3;  //temperature of calculated virial coefficient
		double simTemperature = 1.0*temperature; //temperature governing sampling of configurations
				
		SimulationGraphic sim = new SimulationGraphic(new Space3D());
		final Phase phase = new Phase(sim);
		phase.setBoundary(sim.space.makeBoundary(Space3D.Boundary.NONE));	
		SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
		species.setNMolecules(nMolecules);
		sim.elementCoordinator.go();
		PairSet pairs = new PairSet(((AtomTreeNodeGroup)phase.getAgent(species).node).childList);
		
		P2Cluster p2 = new P2Cluster(sim.hamiltonian.potential, pairs);
		p2.setTemperature(simTemperature);
		P2LennardJones p2LJ = new P2LennardJones(p2);
		MayerModified f = new MayerModified(p2LJ, sigmaHSMod);
		Cluster ring = new etomica.virial.cluster.Ring(nMolecules, 1.0, f);
		p2.setCluster(ring);
		
		Controller controller = new Controller(sim);		
		DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(sim);
		IntegratorMC integrator = new IntegratorMC(sim);
		integrator.setSleepPeriod(1);
//		integrator.setDoSleep(false);
		integrator.setTemperature(simTemperature);
		MyMCMoveAtom mcMoveAtom = new MyMCMoveAtom(integrator);
		
		MeterVirialB3 meterVirial = new MeterVirialB3(sim, pairs, sigmaHSRef, p2, p2LJ);
		meterVirial.setTemperature(temperature);
		
		species.setDiameter(sigmaHSRef);
		
		DisplayPhase display = new DisplayPhase();
		ColorSchemeByType.setColor(species, java.awt.Color.green);
		
		sim.elementCoordinator.go();
		
		Space3D.Vector origin = new Space3D.Vector(5.,5.,5.);
		SpeciesAgent speciesAgent = phase.getAgent(species);
		speciesAgent.coord.translateTo(new Space3D.Vector(5.,5.,5.));
		speciesAgent.firstMolecule().coord.translateTo(new Space3D.Vector(5.,5.,5.));
				
		AtomList childList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
		AtomIteratorList list1 = new AtomIteratorList(childList);
		list1.reset();
		while(list1.hasNext()) list1.next().coord.translateTo(origin);

		MeterDatumSourceWrapper b3Meter = new MeterDatumSourceWrapper(meterVirial);
		b3Meter.setHistorying(true);
		DisplayPlot bPlot = new DisplayPlot(sim);
		bPlot.setDataSources(b3Meter.getHistory());
		bPlot.setWhichValue(MeterAbstract.CURRENT);
		b3Meter.getHistory().setNValues(1000);
		bPlot.setLabel("B3 running average");
		
		DisplayPlot clusterPlot = new DisplayPlot(sim);
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
		
		sim.elementCoordinator.go();

		sim.makeAndDisplayFrame();		
	
	}
}
