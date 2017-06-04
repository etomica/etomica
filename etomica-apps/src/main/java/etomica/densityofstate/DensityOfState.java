/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.densityofstate;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.IAtomType;
import etomica.box.Box;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.data.histogram.HistogramExpanding;
import etomica.yukawa.P2Yukawa;

/**
 * A Yukawa Monte-Carlo simulation in 3D
 */

public class DensityOfState extends Simulation{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public IntegratorMC integrator;
	public MCMoveAtom mcMoveAtom;
	public SpeciesSpheresMono species;
	public Box box;
	public P2Yukawa potential;
	public Controller controller;
	public ActivityIntegrate activityIntegrate; 
	
	public DensityOfState(){
		this(500);
	}
	
	
	
	public DensityOfState(int numAtoms){
		super(Space3D.getInstance());
		
		
		potentialMaster = new PotentialMasterMonatomic(this);
		integrator = new IntegratorMC(this, potentialMaster);
		mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
		mcMoveAtom.setAtomSource(new AtomSourceRandomLeaf());
		mcMoveAtom.setStepSize(0.2);
		integrator.getMoveManager().addMCMove(mcMoveAtom);
		integrator.getMoveManager().setEquilibrating(false);
		activityIntegrate = new ActivityIntegrate(integrator);
		getController().addAction(activityIntegrate);
		species = new SpeciesSpheresMono(this, space);
		box.setNMolecules(species, numAtoms);
		box = new Box(space);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.65);
        inflater.actionPerformed();
		potential = new P2Yukawa(space);
		double truncationRadius = 3.0*potential.getKappa();
		if(truncationRadius > 0.5*box.getBoundary().getBoxSize().getX(0)){
			throw new RuntimeException("Truncaiton radius too large.  Max allowed is "+0.5*box.getBoundary().getBoxSize().getX(0));
		}
		P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
		((PotentialMasterCell)potentialMaster).setCellRange(3);
		((PotentialMasterCell)potentialMaster).setRange(potentialTruncated.getRange());
		potentialMaster.addPotential(potentialTruncated, new IAtomType[] {species.getLeafType(), species.getLeafType()});
			
		integrator.getMoveEventManager().addListener(((PotentialMasterCell)potentialMaster).getNbrCellManager(box).makeMCMoveListener());
		
		new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
		integrator.setBox(box);
		
		((PotentialMasterCell)potentialMaster).getNbrCellManager(box).assignCellAll();
		
	}
	public PotentialMaster potentialMaster;
	
	public static void main(String[] args){
		
		int numAtoms = 500;
		long maxSteps = 10000000;
		double temperature = 6.0;
		
		if (args.length > 0){
			numAtoms = Integer.valueOf(args[0]).intValue();
		}
		if (args.length > 1){
			temperature = Double.parseDouble(args[1]);
		}
		if (args.length > 2) {
			maxSteps = Long.parseLong(args[2]);
		}
		 
		DensityOfState sim = new DensityOfState(numAtoms);
		sim.activityIntegrate.setMaxSteps(maxSteps);
		sim.integrator.setTemperature(temperature);
		sim.getController().actionPerformed();
				
		MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
		AccumulatorHistogram histogram = new AccumulatorHistogram(new HistogramExpanding(1), 100); //bin size =1
		
		histogram.setPushInterval(100);
		DataPump energyManager = new DataPump(energyMeter, histogram);
		

		sim.activityIntegrate.setMaxSteps(maxSteps);
		sim.integrator.getEventManager().addListener(new IntegratorListenerAction(energyManager));
		sim.getController().reset();
		sim.getController().actionPerformed();
		
		DataLogger dataLogger = new DataLogger();
		DataTableWriter dataTableWriter = new DataTableWriter();
		dataLogger.setFileName("histogram" + "@" + (int)(temperature*10) + ".dat");
		dataLogger.setDataSink(dataTableWriter);
		dataLogger.putDataInfo(histogram.getDataInfo());
		
		dataLogger.setWriteInterval(1);
		dataLogger.putData(histogram.getData());
		dataLogger.closeFile();

		
		
		
		
		//DisplayPlot HistogramDisplay = new DisplayPlot();
		//histogram.addDataSink(HistogramDisplay.getDataSet().makeDataSink());
		
		//SimulationGraphic simGraphic = new SimulationGraphic(sim);
        //simGraphic.makeAndDisplayFrame();
        //simGraphic.panel().add(HistogramDisplay.graphic());
        //ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        //colorScheme.setColor(sim.species.getMoleculeType(), java.awt.Color.red);
		
	}
}
