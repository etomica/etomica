/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.chem.elements.Silver;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.IDataInfo;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

/**
 * @author ub2092
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class MEAM_MC extends Simulation {
	
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM MC";
    public final PotentialMaster potentialMaster;
	public IntegratorMC integrator;
	public SpeciesSpheresMono sn;
    public SpeciesSpheresMono ag;
    public SpeciesSpheresMono cu;
	public Box box;
	public PotentialMEAM potentialN;
	public Controller controller;
	public DisplayPlot plot;
	public MeterEnergy energy;
	public ActivityIntegrate activityIntegrate;
	public IDataInfo info2;

	public MEAM_MC() {
		super(Space3D.getInstance()); //INSTANCE); kmb change 8/3/05
		potentialMaster = new PotentialMaster();
		box = new Box(space);
		integrator = new IntegratorMC(this, potentialMaster, box);
		integrator.getMoveManager().addMCMove(new MCMoveAtom(getRandom(), potentialMaster,
				space, 0.1, 0.2, false));
		integrator.setTemperature(Kelvin.UNIT.toSim(298));
		//integrator.setThermostatInterval(10);
		integrator.setIsothermal(true);
		activityIntegrate = new ActivityIntegrate(integrator);
		activityIntegrate.setSleepPeriod(2);
		getController().addAction(activityIntegrate);
		sn = new SpeciesSpheresMono(space, Tin.INSTANCE);
		ag = new SpeciesSpheresMono(space, Silver.INSTANCE);
		cu = new SpeciesSpheresMono(space, Copper.INSTANCE);

		addSpecies(sn);
		addSpecies(ag);
		addSpecies(cu);
		addBox(box);
		box.setNMolecules(sn, 216);
		box.setNMolecules(ag, 0);
		box.setNMolecules(cu, 0);

		//beta-Sn box

		//The dimensions of the simulation box must be proportional to those of
		//the unit cell to prevent distortion of the lattice.  The values for the
		//lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815
		//angstroms) are taken from the ASM Handbook.
		box.getBoundary().setBoxSize(new Vector3D(5.8314 * 3, 5.8314 * 3, 3.1815 * 6));
		PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
		//Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
		//box.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
		//PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
		BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());

		//FCC Cu
		/**
		 box.setDimensions(new Vector3D(3.6148*3, 3.6148*3, 3.6148*6));
		 PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
		 LatticeCrystal crystal = new LatticeCrystal(new Crystal(
		 primitive, new BasisCubicFcc(primitive)));
		 **/

		//FCC Ag
		/**
		 box.setDimensions(new Vector3D(4.0863*3, 4.0863*3, 4.0863*6));
		 PrimitiveCubic primitive = new PrimitiveCubic(space, 4.0863);
		 LatticeCrystal crystal = new LatticeCrystal(new Crystal(
		 primitive, new BasisCubicFcc(primitive)));
		 **/

		Configuration config = new ConfigurationLattice(crystal, space);
		config.initializeCoordinates(box);

		potentialN = new PotentialMEAM(space);
		potentialN.setParameters(sn.getLeafType(), ParameterSetMEAM.Sn);
		potentialN.setParameters(ag.getLeafType(), ParameterSetMEAM.Ag);
		potentialN.setParameters(cu.getLeafType(), ParameterSetMEAM.Cu);
		potentialN.setParametersIMC(cu.getLeafType(), ParameterSetMEAM.Cu3Sn);
		potentialN.setParametersIMC(ag.getLeafType(), ParameterSetMEAM.Ag3Sn);
		this.potentialMaster.addPotential(potentialN, new AtomType[]{sn.getLeafType(), ag.getLeafType(), cu.getLeafType()});
		BoxImposePbc imposepbc = new BoxImposePbc(space);
		imposepbc.setBox(box);
		integrator.getEventManager().addListener(new IntegratorListenerAction(imposepbc));

		// IntegratorCoordConfigWriter - Displacement output (3/1/06 - MS)
		//IntegratorCoordConfigWriter coordWriter = new IntegratorCoordConfigWriter(space, "MEAMoutput");
		//coordWriter.setBox(box);
		//coordWriter.setIntegrator(integrator);
		//coordWriter.setWriteInterval(100);

		// Control simulation lengths
		//activityIntegrate.setMaxSteps(500);

		energy = new MeterEnergy(potentialMaster, box);
	}

    public static void main(String[] args) {
        MEAM_MC sim = new MEAM_MC();
        MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
        energyMeter.setBox(sim.box);
        AccumulatorHistory energyAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        final DisplayPlot plot = new DisplayPlot();
        energyAccumulator.setDataSink(plot.getDataSet().makeDataSink());
        DataPump energyPump = new DataPump(energyMeter, energyAccumulator);
        //energyAccumulator.setBlockSize(50);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(energyPump));

        SimulationGraphic simgraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, APP_NAME);
        simgraphic.getController().getDataStreamPumps().add(energyPump);

        simgraphic.getPanel().plotPanel.add(plot.graphic(), SimulationPanel.getVertGBC());
        //simgraphic.panel().add(plotKE.graphic());

        simgraphic.getController().getReinitButton().setPostAction(simgraphic.getPaintAction(sim.box));

        ColorSchemeByType colorScheme = ((ColorSchemeByType) ((DisplayBox) simgraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.sn.getLeafType(), java.awt.Color.blue);
        colorScheme.setColor(sim.ag.getLeafType(), java.awt.Color.gray);
        colorScheme.setColor(sim.cu.getLeafType(), java.awt.Color.orange);

        simgraphic.makeAndDisplayFrame(APP_NAME);

        //sim.activityIntegrate.setMaxSteps(1000);
        //sim.getController().run();
        //DataGroup data = (DataGroup)energyAccumulator.getData(); // kmb change type to Data instead of double[]
        //double PE = ((DataDouble)data.getData(AccumulatorAverage.StatType.AVERAGE.index)).x
        // /sim.species.getAgent(sim.box).getNMolecules();  // kmb changed 8/3/05
        //double PE = data[AccumulatorAverage.AVERAGE.index]/sim.species.getAgent(sim.box).getNMolecules();  // orig line
        //System.out.println("PE(eV)="+ElectronVolt.UNIT.fromSim(PE));
    }

}
