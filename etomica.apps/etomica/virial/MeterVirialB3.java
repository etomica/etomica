package etomica.virial;

import etomica.*;
import etomica.virial.cluster.C3;
import etomica.units.Dimension;
import etomica.graphics.*;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class MeterVirialB3 extends MeterGroup implements DatumSource {

	private Atom[] atoms;
	public static final int N = 3;
	private double sigma, sigma2;
	private double[][] f, fSim, fRef;
	private AtomPair[][] pairs;
	private Cluster[] clusters;
	private double B3Ref;
	private P2MayerModified simulationPotential;
	private double temperature, beta;

	/**
	 * Constructor for MeterVirialB3.
	 * @param sim
	 * @param nMeters
	 */
	public MeterVirialB3(Simulation sim) {
		super(sim, 2);
		atoms = new Atom[N];
		f = new double[N-1][];
		fSim = new double[N-1][];
		fRef = new double[N-1][];
		pairs = new AtomPair[N-1][];
		for(int i=0; i<N-1; i++) {
			f[i] = new double[N-1-i];
			fSim[i] = new double[N-1-i];
			fRef[i] = new double[N-1-i];
			pairs[i] = new AtomPair[N-1-i];
			for(int j=0; j<N-1-i; j++) pairs[i][j] = new AtomPair(sim.space);
		}
		clusters = new Cluster[1];
		clusters[0] = new C3();
		setSigma(Default.ATOM_SIZE);
	}

	public void setAtoms(AtomList list) {
		AtomIteratorList iterator = new AtomIteratorList(list);
		iterator.reset();
		int i=0;
		while(iterator.hasNext() && i<N) atoms[i++] = iterator.next();
		if(atoms[N-1]==null) throw new RuntimeException("Not enough atoms given to MeterVirial");
		for(i=0; i<N-1; i++) {
			for(int j=0; j<N-1-i; j++) pairs[i][j].reset(atoms[i],atoms[i+j+1]);
		}
	}
	

	/**
	 * @see etomica.MeterGroup#updateValues()
	 */
	public void updateValues() {
		for(int i=0; i<N-1; i++) {
			for(int j=0; j<N-1-i; j++) {
				pairs[i][j].reset();
				fSim[i][j] = Math.exp(-beta*simulationPotential.energy(pairs[i][j]));
				double bu = simulationPotential.mostRecentBetaU();
				f[i][j] = Math.exp(-bu) - 1.0;
				fRef[i][j] = (pairs[i][j].r2()<sigma2) ? -1.0 : 0.0;
			}
		}
		double CSim = clusters[0].value(fSim);  //ring
		double C3 = clusters[0].value(f);
		double C3Ref = clusters[0].value(fRef);
		currentValues[1] = C3Ref/CSim;
		currentValues[0] = C3/CSim;
	}

	/**
	 * @see etomica.MeterAbstract#getDimension()
	 */
	public Dimension getDimension() {
		return null;
	}

	/**
	 * @see etomica.DatumSource#value(etomica.DataSource.ValueType)
	 */
	public double value(DataSource.ValueType type) {
		MeterScalar[] m = allMeters();
		double sum = 0.0;
		for(int i=0; i<clusters.length; i++) {
			sum += clusters[i].weight()*m[i].average();
		}
		double denom = clusters[0].weight() * m[1].average();
		return (denom != 0.0) ? B3Ref*sum/denom : Double.NaN;
	}

	public static void main(String[] args) {
	
		Default.makeLJDefaults();
		Default.TRUNCATE_POTENTIALS = false;
//		Default.EXPLICIT_LOOP = true;
		
		SimulationGraphic sim = new SimulationGraphic(new Space3D());
		final Phase phase = new Phase(sim);
		phase.setBoundary(sim.space.makeBoundary(Space3D.Boundary.NONE));	
		SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
		species.setNMolecules(3);
		Controller controller = new Controller(sim);		
		DeviceTrioControllerButton controlPanel = new DeviceTrioControllerButton(sim);
		IntegratorMC integrator = new IntegratorMC(sim);
		integrator.setSleepPeriod(1);
		integrator.setDoSleep(false);
		MyMCMoveAtom mcMoveAtom = new MyMCMoveAtom(integrator);
		
		P2MayerModified p2Mayer = new P2MayerModified(sim.hamiltonian.potential);
		P2LennardJones p2LJ = new P2LennardJones(p2Mayer);
		p2Mayer.setPotential(p2LJ);
		p2Mayer.setSpecies(species, species);
		
		MeterVirialB3 meterVirial = new MeterVirialB3(sim);
		
		double sigmaHS = 1.28412293285;  //  ( 4 + 4 sqrt(1-Ln(2)) ) / Ln(4))^(1/6), which is where f(r) = 1 for LJ
		double temperature = 4.3;
		
		p2Mayer.setSigma(sigmaHS);
		species.setDiameter(sigmaHS);
		
		p2Mayer.setTemperature(temperature);
		integrator.setTemperature(temperature);
		
		DisplayPhase display = new DisplayPhase();
		ColorSchemeByType.setColor(species, java.awt.Color.green);
		

		AtomIteratorNeighbor nbrIterator1 = new AtomIteratorNeighbor(true);
		sim.elementCoordinator.go();
		
		Space3D.Vector origin = new Space3D.Vector(5.,5.,5.);
		SpeciesAgent speciesAgent = phase.getAgent(species);
		speciesAgent.coord.translateTo(new Space3D.Vector(5.,5.,5.));
		speciesAgent.firstMolecule().coord.translateTo(new Space3D.Vector(5.,5.,5.));
		
		meterVirial.setSigma(sigmaHS);
		meterVirial.setAtoms(((AtomTreeNodeGroup)speciesAgent.node).childList);
		meterVirial.setTemperature(temperature);
		meterVirial.setSimulationPotential(p2Mayer);
		
		NeighborManager.Criterion criterion = new NeighborManager.Criterion() {
			public boolean areNeighbors(Atom a1, Atom a2) {
				return Math.abs(a1.node.index()-a2.node.index()) == 1 ||
				   (a1==phase.firstAtom() && a2==phase.lastAtom()) ||
				   (a2==phase.firstAtom() && a1==phase.lastAtom());
			}};
//				return Math.abs(a1.node.index()-a2.node.index()) == 1;
//			}};
		
		AtomList childList = ((AtomTreeNodeGroup)phase.getAgent(species).node).childList;
		nbrIterator1.setupNeighbors(childList, criterion);
		AtomIteratorList list1 = new AtomIteratorList(childList);
		AtomPairIterator api1= new ApiGeneral(sim.space,list1,nbrIterator1);
		p2Mayer.setIterator(api1);
		list1.reset();
		while(list1.hasNext()) list1.next().coord.translateTo(origin);
//		IteratorDirective.testSuitePair(api1,phase.firstAtom(),phase.firstAtom(),phase.lastAtom());

		MeterDatumSourceWrapper b4Meter = new MeterDatumSourceWrapper(meterVirial);
		b4Meter.setHistorying(true);
		DisplayPlot bPlot = new DisplayPlot(sim);
		bPlot.setDataSources(b4Meter.getHistory());
		bPlot.setWhichValue(MeterAbstract.CURRENT);
		b4Meter.getHistory().setNValues(1000);
		bPlot.setLabel("B3 running average");
		
		DisplayPlot clusterPlot = new DisplayPlot(sim);
		meterVirial.setHistorying(true);
		MeterDatumSourceWrapper[] clusterMeter = new MeterDatumSourceWrapper[4];
		for(int i=0; i<2; i++) {
			clusterMeter[i] = new MeterDatumSourceWrapper(meterVirial.allMeters()[i]);
			clusterMeter[i].setHistorying(true);
			clusterMeter[i].getHistory().setNValues(1000);
		}
		clusterMeter[0].setLabel("Reference");
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
	/**
	 * Returns the sigma.
	 * @return double
	 */
	public double getSigma() {
		return sigma;
	}

	/**
	 * Sets the sigma.
	 * @param sigma The sigma to set
	 */
	public void setSigma(double sigma) {
		this.sigma = sigma;
		sigma2 = sigma*sigma;
		double b0 = 2.0*Math.PI/3.0*sigma2*sigma;
		B3Ref = 5./8.*b0*b0;
	}
	
	/**
	 *  Overrides MCMoveAtom to prevent index-0 molecule from being displaced
	 */
	private static class MyMCMoveAtom extends MCMoveAtom {
		MyMCMoveAtom(IntegratorMC integrator) {super(integrator);}
		
		public boolean doTrial() {
			if(phase.atomCount()==0) return false;
			atom = null;
			while(atom == null || atom.node.index()==0) atom = phase.speciesMaster.atomList.getRandom();
			uOld = potential.calculate(phase, iteratorDirective.set(atom), energy.reset()).sum();
			atom.coord.displaceWithin(stepSize);
			uNew = Double.NaN;
			return true;
		}
	}

	/**
	 * Returns the temperature.
	 * @return double
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * Sets the temperature.
	 * @param temperature The temperature to set
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
		beta = 1.0/temperature;
	}

	/**
	 * Returns the simulationPotential.
	 * @return P2MayerModified
	 */
	public P2MayerModified getSimulationPotential() {
		return simulationPotential;
	}

	/**
	 * Sets the simulationPotential.
	 * @param simulationPotential The simulationPotential to set
	 */
	public void setSimulationPotential(P2MayerModified simulationPotential) {
		this.simulationPotential = simulationPotential;
	}


}
