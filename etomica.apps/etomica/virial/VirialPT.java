package etomica.virial;

import etomica.*;
import etomica.graphics.*;

public class VirialPT extends SimulationGraphic {
    
	Phase[] phase;
	SpeciesSpheresMono species;
	Controller controller;
	IntegratorMC[] integrator;
	P2Cluster potential;
    
    
	public VirialPT(int nM, double betaMin, double betaMax, int nPhase) {
		super(new Space3D());
		Simulation.instance = this;
		Default.TRUNCATE_POTENTIALS = false;
		Default.TEMPERATURE = 1.0; 
		Default.HISTORY_PERIOD = 1000; //for phase tracking
        
		elementCoordinator.addMediatorPair(new MediatorGraphic.DisplayMeter.NoAction(elementCoordinator));
        
		species = new SpeciesSpheresMono();
		species.setNMolecules(nM);
		
		P2LennardJones p2LJ = new P2LennardJones(Simulation.instance.hamiltonian.potential);//parent of this potential will not be connected to the simulation
		MayerGeneral f = new MayerGeneral(p2LJ);
		
		//take care to define other clusters (below) appropriately if using ReeHoover
//		Cluster targetCluster = new Full(nMolecules, 1.0, f);
		Cluster targetCluster = new Cluster(5,1.0,
			new Cluster.BondGroup(f, new int[][] {{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{2,3},{2,4},{3,4}}));
//		Cluster targetCluster = new etomica.virial.cluster.ReeHoover(4, 1.0, new Cluster.BondGroup(f, Standard.D4));
		P2ClusterSigned potential = new P2ClusterSigned(refSimulation.hamiltonian.potential, pairs);

		AtomType type = new AtomType.SphereVariable(species.moleculeFactory(), potential.index);
		((AtomFactoryMono)species.moleculeFactory()).setType(type);
		species.protoType = (AtomType.Sphere)type;
        
		IntegratorPT.MCMoveSwapFactory swapFactory = new MCMoveSwapPolydisperse.Factory(potential.index);
		IntegratorPT integratorPT = new IntegratorPT(this, swapFactory);
		integratorPT.setDoSleep(false);
		integratorPT.setInterval(species.getNMolecules());
		integratorPT.setSwapInterval(1000);
        
		InfPoly muFunction = new InfPoly(c0Min);
		phase = new Phase[nPhase];
		integrator = new IntegratorMC[nPhase];
		MeterAcceptanceRatio[] meterAccept = new MeterAcceptanceRatio[nPhase-1];
		MeterM0[] meterM0 = new MeterM0[nPhase];
		MeterPerturbToPure[] meterPerturb = new MeterPerturbToPure[nPhase];
		double c0Last = 0.0;
		for(int i=0; i<nPhase; i++) {
			double c0; 
			if(nPhase > 1) c0 = -1 + Math.exp(Math.log(c0Min+1) + (double)i/(double)(nPhase-1)*Math.log((c0Max+1)/(c0Min+1)));
			else c0 = c0Min;
            
			phase[i] = new Phase();
			integrator[i] = new IntegratorMC(space, hamiltonian.potential);
			phase[i].setIntegrator(integrator[i]);
			integrator[i].addIntervalListener(new PhaseAction.ImposePbc(phase[i]));
            
			MCMoveAtom moveAtom = new MCMoveAtom(integrator[i]);
			MCMoveDiameter moveDiameter = new MCMoveDiameter(integrator[i], potential.index);
			muFunction = new InfPoly(c0);
			moveDiameter.setMu(muFunction);
            
			integratorPT.addIntegrator(integrator[i]);
 
			DisplayPhase display = new DisplayPhase();
			display.setPhase(phase[i]);
			display.setLabel(Double.toString(c0));
            
			if(i>0) {
				meterAccept[i-1] = new MeterAcceptanceRatio(this, integratorPT.swapMoves()[i-1]);
				meterAccept[i-1].setLabel(Double.toString(c0Last)+"<-->"+Double.toString(c0));
			}
			c0Last = c0;
            
			meterM0[i] = new MeterM0(this, potential.index);
			meterM0[i].setPhase(phase[i]);
			meterM0[i].setLabel(Double.toString(c0));
            
			meterPerturb[i] = new MeterPerturbToPure(this, species, potential.index, c0);
			meterPerturb[i].setPhase(phase[i]);
			meterPerturb[i].setLabel(Double.toString(c0));
		}
        
		IntegratorPT.MeterPhaseTracker phaseTracker = integratorPT.new MeterPhaseTracker();
		phaseTracker.setPhases(phase);
		phaseTracker.setUpdateInterval(10);
		DisplayPlot phaseTrackPlot = new DisplayPlot();
		phaseTrackPlot.setDoLegend(false);
		phaseTrackPlot.setWhichValue(MeterAbstract.CURRENT);
		phaseTrackPlot.setDataSources(phaseTracker.dataSource());
		phaseTrackPlot.setLabel("Track");
		phaseTrackPlot.setUpdateInterval(10);
        
		DisplayTable acceptanceTable = new DisplayTable();
		acceptanceTable.setDatumSources(meterAccept);
		acceptanceTable.setLabel("Accept");
		acceptanceTable.setWhichValues(MeterAbstract.CURRENT);
        
		DisplayTable m0Table = new DisplayTable();
		m0Table.setDatumSources(meterM0);
		m0Table.setLabel("m0");
		m0Table.setWhichValues(new DataSource.ValueType[] 
						{MeterAbstract.CURRENT, MeterAbstract.AVERAGE, MeterAbstract.ERROR});

		DisplayTableFunction perturbTable = new DisplayTableFunction();
		perturbTable.setPrecision(7);
		perturbTable.setDataSources(meterPerturb);
		perturbTable.setLabel("Perturb");
		perturbTable.setWhichValue(MeterAbstract.AVERAGE);
		DisplayTableFunction perturbTableError = new DisplayTableFunction();
		perturbTableError.setDataSources(meterPerturb);
		perturbTableError.setLabel("Perturb Err");
		perturbTableError.setWhichValue(MeterAbstract.ERROR);
        
		DeviceButton logTableButton = new DeviceButton(perturbTable.makeLogTableAction());
		logTableButton.setLabel("Perturb table out");
        
		DeviceButton logTableButton1 = new DeviceButton(m0Table.makeLogTableAction());
		logTableButton1.setLabel("m0 table out");
        
                        
		controller = new Controller();
        
		DisplayBox cyclesBox = new DisplayBox(new MeterCycles());
		cyclesBox.setPrecision(6);
        
		ColorSchemeByType.setColor(species,java.awt.Color.green);
		DeviceTrioControllerButton buttons = new DeviceTrioControllerButton();
		if(nPhase == 1) {
			DeviceSlider slider = new DeviceSlider(new Modulator(muFunction, "c0"));
			slider.setMinimum(-1);
			slider.setMaximum(29);
		}
	}
    
	public class InfPoly implements etomica.utility.Function {
		private double c0;
		public InfPoly(double c0) {this.c0 = c0;}
        
		public double f(double x) {return c0*Math.log(x);}
		public double inverse(double f) {return Math.exp(f/c0);}
		public double dfdx(double x) {return c0/x;}
		public void setC0(double c0) {this.c0 = c0;}
		public double getC0() {return c0;}
	}
    
	public static void main(String[] args) {
		int nM = Integer.parseInt(args[0]);
		double c0Min = Double.parseDouble(args[1]);
		double c0Max = Double.parseDouble(args[2]);
		int n = Integer.parseInt(args[3]);
		PolydisperseHS sim = new PolydisperseHS(nM, c0Min, c0Max, n);
		sim.elementCoordinator.go();
		sim.makeAndDisplayFrame();
	}
}