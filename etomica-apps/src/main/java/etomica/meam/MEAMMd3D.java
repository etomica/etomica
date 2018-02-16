/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.chem.elements.Silver;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.*;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.CriterionSimple;
import etomica.nbr.list.PotentialMasterList;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.*;

import java.util.ArrayList;

/**
 * Molecular-Dynamics Simulation Using the Modified Embedded-Atom Method 
 * (MEAM) Potential.  
 * 
 * The MEAM potential is intended for use with metallic and covalently-bonded
 * solid systems.
 * 
 * The MEAM potential for an atom is built using terms describing parts of the
 * relationships between the atom and each of its neighbors, the number of which 
 * is determined by a cutoff and/or screening function.  Each type of pair-
 * wise term is summed over all the neighbors, and then used in expressions 
 * describing the embedding energy and the repulsive energy of the atom.   
 * Effectively, the MEAM potential is a many-body potential.  
 * 
 * This class was adapted from LjMd3D.java by K.R. Schadel and A. Schultz in July 
 * 2005.  Intitially, it employed a version of the embedded-atom method potential, 
 * and was later adapted in February 2006 to use the modified embedded-atom method
 * potential.
 */
 
public class MEAMMd3D extends Simulation {
    
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM Md3D";
    public final PotentialMasterList potentialMaster;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono sn;
    public SpeciesSpheresMono ag;
    public SpeciesSpheresMono cu;
    public Box box;
    public PotentialMEAM potentialN;
    public Controller controller;
    public DisplayBox display;
    public DisplayPlot plot;
    public MeterEnergy energy;
    public ActivityIntegrate activityIntegrate;
    public IDataInfo info2;

    public MEAMMd3D() {
        super(Space3D.getInstance()); //INSTANCE); kmb change 8/3/05
        potentialMaster = new PotentialMasterList(this, space);
        box = new Box(space);
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setTimeStep(0.001);
        integrator.setTemperature(Kelvin.UNIT.toSim(295));
        integrator.setThermostatInterval(100);
        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        sn = new SpeciesSpheresMono(space, Tin.INSTANCE);
        sn.setIsDynamic(true);
        ag = new SpeciesSpheresMono(space, Silver.INSTANCE);
        ag.setIsDynamic(true);
        cu = new SpeciesSpheresMono(space, Copper.INSTANCE);
        cu.setIsDynamic(true);

        addSpecies(sn);
        addSpecies(ag);
        addSpecies(cu);
        addBox(box);
        box.setNMolecules(sn, 0);
        box.setNMolecules(ag, 256);
        box.setNMolecules(cu, 0);

        // beta-Sn box
        /**
         //The dimensions of the simulation box must be proportional to those of
         //the unit cell to prevent distortion of the lattice.  The values for the
         //lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815
         //angstroms) are taken from the ASM Handbook.
         box.setDimensions(new Vector3D(5.8314*3, 5.8314*3, 3.1815*6));
         PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
         //Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
         //box.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
         //PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
         LatticeCrystal crystal = new LatticeCrystal(new Crystal(
         primitive, new BasisBetaSnA5(primitive)));
         */

        //FCC Cu
        /**
         box.setDimensions(new Vector3D(3.6148*4, 3.6148*4, 3.6148*4));
         PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
         LatticeCrystal crystal = new LatticeCrystal(new Crystal(
         primitive, new BasisCubicFcc(primitive)));
         */

        //FCC Ag

        box.getBoundary().setBoxSize(new Vector3D(4.0863 * 4, 4.0863 * 4, 4.0863 * 4));
        PrimitiveCubic primitive = new PrimitiveCubic(space, 4.0863);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());


        Configuration config = new ConfigurationLattice(crystal, space);
        config.initializeCoordinates(box);

        potentialN = new PotentialMEAM(space);
        potentialN.setParameters(sn.getLeafType(), ParameterSetMEAM.Sn);
        potentialN.setParameters(ag.getLeafType(), ParameterSetMEAM.Ag);
        potentialN.setParameters(cu.getLeafType(), ParameterSetMEAM.Cu);
        potentialN.setParametersIMC(cu.getLeafType(), ParameterSetMEAM.Cu3Sn);
        potentialN.setParametersIMC(ag.getLeafType(), ParameterSetMEAM.Ag3Sn);
        this.potentialMaster.addPotential(potentialN, new AtomType[]{sn.getLeafType(), ag.getLeafType(), cu.getLeafType()});
        potentialMaster.setRange(potentialN.getRange() * 1.1);
        potentialMaster.setCriterion(potentialN, new CriterionSimple(this, space, potentialN.getRange(), potentialN.getRange() * 1.1));
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

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
        MEAMMd3D sim = new MEAMMd3D();

        MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
        MeterKineticEnergy kineticMeter = new MeterKineticEnergy(sim.box);

        energyMeter.setBox(sim.box);

        AccumulatorHistory energyAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorHistory kineticAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());

        AccumulatorAverageCollapsing accumulatorAveragePE = new AccumulatorAverageCollapsing();
        AccumulatorAverageCollapsing accumulatorAverageKE = new AccumulatorAverageCollapsing();

        DataPump energyPump = new DataPump(energyMeter, accumulatorAveragePE);
        DataPump kineticPump = new DataPump(kineticMeter, accumulatorAverageKE);

        accumulatorAveragePE.addDataSink(energyAccumulator, new StatType[]{AccumulatorAverage.MOST_RECENT});
        accumulatorAverageKE.addDataSink(kineticAccumulator, new StatType[]{AccumulatorAverage.MOST_RECENT});

        DisplayPlot plotPE = new DisplayPlot();
        plotPE.setLabel("PE Plot");
        DisplayPlot plotKE = new DisplayPlot();
        plotKE.setLabel("KE Plot");

        energyAccumulator.setDataSink(plotPE.getDataSet().makeDataSink());
        kineticAccumulator.setDataSink(plotKE.getDataSet().makeDataSink());
        //energyAccumulator.setBlockSize(50);

        accumulatorAveragePE.setPushInterval(1);
        accumulatorAverageKE.setPushInterval(1);

        //Heat Capacity (PE)
        DataProcessorCvMD dataProcessorPE = new DataProcessorCvMD();
        dataProcessorPE.setIntegrator(sim.integrator);

        //Heat Capacity (KE)
        DataProcessorCvMD dataProcessorKE = new DataProcessorCvMD();
        dataProcessorKE.setIntegrator(sim.integrator);

        accumulatorAveragePE.addDataSink(dataProcessorPE, new StatType[]{AccumulatorAverage.STANDARD_DEVIATION});
        accumulatorAverageKE.addDataSink(dataProcessorKE, new StatType[]{AccumulatorAverage.STANDARD_DEVIATION});

        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(energyPump));
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(kineticPump));

        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
        ArrayList dataStreamPumps = simGraphic.getController().getDataStreamPumps();
        dataStreamPumps.add(energyPump);
        dataStreamPumps.add(kineticPump);

        DisplayTextBox cvBoxPE = new DisplayTextBox();
        dataProcessorPE.setDataSink(cvBoxPE);
        cvBoxPE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double[]{1, -1, -1}));
        cvBoxPE.setLabel("PE Cv contrib.");
        DisplayTextBox cvBoxKE = new DisplayTextBox();
        dataProcessorKE.setDataSink(cvBoxKE);
        cvBoxKE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double[]{1, -1, -1}));
        cvBoxKE.setLabel("KE Cv contrib.");

        simGraphic.add(/*"PE Plot",*/plotPE);
        simGraphic.add(/*"KE Plot",*/plotKE);

        simGraphic.getPanel().controlPanel.add(cvBoxKE.graphic(), SimulationPanel.getVertGBC());
        simGraphic.getPanel().controlPanel.add(cvBoxPE.graphic(), SimulationPanel.getVertGBC());

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        ColorSchemeByType colorScheme = ((ColorSchemeByType) ((DisplayBox) simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.sn.getLeafType(), java.awt.Color.blue);
        colorScheme.setColor(sim.ag.getLeafType(), java.awt.Color.gray);
        colorScheme.setColor(sim.cu.getLeafType(), java.awt.Color.orange);

        simGraphic.makeAndDisplayFrame(APP_NAME);

        //sim.activityIntegrate.setMaxSteps(1000);
        //sim.getController().run();
        //DataGroup data = (DataGroup)energyAccumulator.getData(); // kmb change type to Data instead of double[]
        //double PE = ((DataDouble)data.getData(AccumulatorAverage.StatType.AVERAGE.index)).x
        // /sim.species.getAgent(sim.box).getNMolecules();  // kmb changed 8/3/05
        //double PE = data[AccumulatorAverage.AVERAGE.index]/sim.species.getAgent(sim.box).getNMolecules();  // orig line
        //System.out.println("PE(eV)="+ElectronVolt.UNIT.fromSim(PE));
    }
    
}
