/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.chem.elements.Tin;
import etomica.config.GrainBoundaryConfiguration;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.*;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.CriterionSimple;
import etomica.nbr.list.PotentialMasterList;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
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
 
public class MEAM_3DMDwithSnCuGB extends Simulation {
    
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM 3DMD w/SnCuGB";
    public final PotentialMasterList potentialMaster;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono snFixedA;
    public SpeciesSpheresMono snA;
//    public SpeciesSpheresMono agA;
//    public SpeciesSpheresMono cuA;
    public SpeciesSpheresMono cuFixedB;
//    public SpeciesSpheresMono snB;
//    public SpeciesSpheresMono agB;
    public SpeciesSpheresMono cuB;
    public Box box;
    public PotentialMEAM potentialN;
    public Controller controller;
    public DisplayBox display;
    public DisplayPlot plot;
    public MeterEnergy energy;
    public ActivityIntegrate activityIntegrate;
    public IDataInfo info2;

    public MEAM_3DMDwithSnCuGB() {
        super(Space3D.getInstance());//INSTANCE); kmb change 8/3/05
        box = new Box(new BoundaryRectangularSlit(2, space), space);
        addBox(box);
        potentialMaster = new PotentialMasterList(this, space);
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, space, box);
        integrator.setTimeStep(0.001);
        integrator.setTemperature(Kelvin.UNIT.toSim(295));
        integrator.setThermostatInterval(100);
        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        Tin SnF = new Tin("SnF", Double.POSITIVE_INFINITY);
        snFixedA = new SpeciesSpheresMono(space, SnF);
        snFixedA.setIsDynamic(true);
        snA = new SpeciesSpheresMono(space, Tin.INSTANCE);
        snA.setIsDynamic(true);
//        agA = new SpeciesSpheresMono(space, Silver.INSTANCE);
//        cuA = new SpeciesSpheresMono(space, Copper.INSTANCE);
        Copper CuF = new Copper("CuF", Double.POSITIVE_INFINITY);
        cuFixedB = new SpeciesSpheresMono(space, CuF);
        cuFixedB.setIsDynamic(true);
//        snB = new SpeciesSpheresMono(space, Tin.INSTANCE);
//        agB = new SpeciesSpheresMono(space, Silver.INSTANCE);
        cuB = new SpeciesSpheresMono(space, Copper.INSTANCE);
        cuB.setIsDynamic(true);

        addSpecies(snFixedA);
        addSpecies(snA);
//        addSpecies(agA);
//        addSpecies(cuA);
        addSpecies(cuFixedB);
//        addSpecies(snB);
//        addSpecies(agB);
        addSpecies(cuB);


        double aA, bA, cA, aB, bB, cB;
        int nCellsAx, nCellsAy, nCellsAz, nAMobile, nAFixed, basisA,
                nCellsBx, nCellsBy, nCellsBz, nBMobile, nBFixed, basisB,
                nA, nB, nAImpurity, nBImpurity, nAVacancy, nBVacancy;

        nCellsAx = 5;
        nCellsAy = 5;
        nCellsAz = 6;
        nCellsBx = 8;
        nCellsBy = 8;
        nCellsBz = 6;
        nAImpurity = 0;
        nAVacancy = 0;
        nBImpurity = 0;
        nBVacancy = 0;

        // beta-Sn box
        //The values for the lattice parameters for tin's beta box
        //(a = 5.8314 angstroms, c = 3.1815
        //angstroms) are taken from the ASM Handbook.
        aA = bA = 5.8;
        cA = 3.1815;
        basisA = 4;
        PrimitiveTetragonal primitiveA = new PrimitiveTetragonal(space, aA, cA);
        //Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
        //box.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
        //PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
        BravaisLatticeCrystal latticeA = new BravaisLatticeCrystal(primitiveA, new BasisBetaSnA5());
        //FCC Cu
        aB = bB = cB = 3.625;
        basisB = 4;
        PrimitiveCubic primitiveB = new PrimitiveCubic(space, aB);
        BravaisLatticeCrystal latticeB = new BravaisLatticeCrystal(primitiveB, new BasisCubicFcc());

        box.getBoundary().setBoxSize(new Vector3D(aA * nCellsAx, aA * nCellsAy, (cA * nCellsAz) + (cB * nCellsBz)));

        nA = (nCellsAx * nCellsAy * nCellsAz) * basisA;
        nAFixed = (nCellsAx * nCellsAy * 2) * basisA;
        nAMobile = nA - nAFixed - nAImpurity - nAVacancy;
        nB = (nCellsBx * nCellsBy * nCellsBz) * basisB;
        nBFixed = (nCellsBx * nCellsBy * 2) * basisB;
        nBMobile = nB - nBFixed - nBImpurity - nBVacancy;

        box.setNMolecules(snFixedA, nAFixed);
        box.setNMolecules(snA, nAMobile);
//	    agA.getAgent(box).setNMolecules(0);
//	    cuA.getAgent(box).setNMolecules(0);
        box.setNMolecules(cuFixedB, nBFixed);
//	    snB.getAgent(box).setNMolecules(0);
//	    agB.getAgent(box).setNMolecules(0);
        box.setNMolecules(cuB, nBMobile);


        GrainBoundaryConfiguration config = new GrainBoundaryConfiguration(latticeA, latticeB, space);
        config.setDimensions(nCellsAx, nCellsAy, nCellsAz, nCellsBx, nCellsBy,
                nCellsBz, aA, bA, cA, aB, bB, cB);
        config.initializeCoordinates(box);

        potentialN = new PotentialMEAM(space);
        potentialN.setParameters(snFixedA.getLeafType(), ParameterSetMEAM.Sn);
        potentialN.setParameters(snA.getLeafType(), ParameterSetMEAM.Sn);
//		potentialN.setParameters(agA.getLeafType(), ParameterSetMEAM.Ag);
//		potentialN.setParameters(cuA.getLeafType(), ParameterSetMEAM.Cu);
        potentialN.setParameters(cuFixedB.getLeafType(), ParameterSetMEAM.Cu);
//		potentialN.setParameters(snB.getLeafType(), ParameterSetMEAM.Sn);
//		potentialN.setParameters(agB.getLeafType(), ParameterSetMEAM.Ag);
        potentialN.setParameters(cuB.getLeafType(), ParameterSetMEAM.Cu);
//		potentialN.setParametersIMC(cuA.getLeafType(), ParameterSetMEAM.Cu3Sn);
//		potentialN.setParametersIMC(agA.getLeafType(), ParameterSetMEAM.Ag3Sn);
        potentialN.setParametersIMC(cuB.getLeafType(), ParameterSetMEAM.Cu3Sn);
        potentialN.setParametersIMC(cuFixedB.getLeafType(), ParameterSetMEAM.Cu3Sn);
//		potentialN.setParametersIMC(agB.getLeafType(), ParameterSetMEAM.Ag3Sn);
        potentialMaster.addPotential(potentialN, new AtomType[]{snFixedA.getLeafType(), snA.getLeafType(), cuFixedB.getLeafType(), cuB.getLeafType()});
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
        MEAM_3DMDwithSnCuGB sim = new MEAM_3DMDwithSnCuGB();

        MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
        MeterKineticEnergy kineticMeter = new MeterKineticEnergy();

        energyMeter.setBox(sim.box);
        kineticMeter.setBox(sim.box);

        AccumulatorHistory energyAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorHistory kineticAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());

        AccumulatorAverageCollapsing accumulatorAveragePE = new AccumulatorAverageCollapsing();
        AccumulatorAverageCollapsing accumulatorAverageKE = new AccumulatorAverageCollapsing();

        DataPump energyPump = new DataPump(energyMeter, accumulatorAveragePE);
        DataPump kineticPump = new DataPump(kineticMeter, accumulatorAverageKE);

        accumulatorAveragePE.addDataSink(energyAccumulator, new StatType[]{AccumulatorAverage.MOST_RECENT});
        accumulatorAverageKE.addDataSink(kineticAccumulator, new StatType[]{AccumulatorAverage.MOST_RECENT});

        DisplayPlot plot = new DisplayPlot();
        //DisplayPlot plotKE = new DisplayPlot();

        energyAccumulator.setDataSink(plot.getDataSet().makeDataSink());
        kineticAccumulator.setDataSink(plot.getDataSet().makeDataSink());
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

        SimulationGraphic simgraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, APP_NAME);
        ArrayList dataStreamPumps = simgraphic.getController().getDataStreamPumps();
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

        simgraphic.getPanel().plotPanel.add(plot.graphic(), SimulationPanel.getVertGBC());
        //simgraphic.panel().add(plotKE.graphic());
        simgraphic.getPanel().plotPanel.add(cvBoxPE.graphic(), SimulationPanel.getVertGBC());
        simgraphic.getPanel().plotPanel.add(cvBoxKE.graphic(), SimulationPanel.getVertGBC());

        simgraphic.getController().getReinitButton().setPostAction(simgraphic.getPaintAction(sim.box));

        ColorSchemeByType colorScheme = ((ColorSchemeByType) ((DisplayBox) simgraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.snFixedA.getLeafType(), java.awt.Color.white);
        colorScheme.setColor(sim.snA.getLeafType(), java.awt.Color.white);
//    	colorScheme.setColor(sim.agA.getMoleculeType(),java.awt.Color.gray);
//    	colorScheme.setColor(sim.cuA.getMoleculeType(),java.awt.Color.orange);
        colorScheme.setColor(sim.cuFixedB.getLeafType(), java.awt.Color.orange);
//    	colorScheme.setColor(sim.snB.getMoleculeType(),java.awt.Color.white);
//    	colorScheme.setColor(sim.agB.getMoleculeType(),java.awt.Color.gray);
        colorScheme.setColor(sim.cuB.getLeafType(), java.awt.Color.orange);

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
