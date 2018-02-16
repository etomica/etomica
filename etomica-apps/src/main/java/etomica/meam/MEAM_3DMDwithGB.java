/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
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

import java.awt.*;
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
 
public class MEAM_3DMDwithGB extends Simulation {
    
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM 3DMD w/GB";
    public PotentialMasterList potentialMaster;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono snFixedA;
    public SpeciesSpheresMono snA;
//    public SpeciesSpheresMono agA;
//    public SpeciesSpheresMono cuA;
    public SpeciesSpheresMono snFixedB;
    public SpeciesSpheresMono snB;
//    public SpeciesSpheresMono agB;
//    public SpeciesSpheresMono cuB;
    public Box box;
    public PotentialMEAM potentialN;
    public Controller controller;
    public DisplayPlot plot;
    public MeterEnergy energy;
    public ActivityIntegrate activityIntegrate;
    public IDataInfo info2;

    public MEAM_3DMDwithGB() {
        super(Space3D.getInstance()); //INSTANCE); kmb change 8/3/05
        potentialMaster = new PotentialMasterList(this, space);
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
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
        snFixedB = new SpeciesSpheresMono(space, SnF);
        snFixedB.setIsDynamic(true);
        snB = new SpeciesSpheresMono(space, Tin.INSTANCE);
        snB.setIsDynamic(true);
//        agB = new SpeciesSpheresMono(space, Silver.INSTANCE);
//        cuB = new SpeciesSpheresMono(space, Copper.INSTANCE);


        double aA, bA, cA, aB, bB, cB;

        int nCellsAx, nCellsAy, nCellsAz, nAMobile, nAFixed,
                nCellsBx, nCellsBy, nCellsBz, nBMobile, nBFixed,
                nA, nB, nAImpurity, nBImpurity, nAVacancy, nBVacancy;

        nCellsAx = 6;
        nCellsAy = 6;
        nCellsAz = 4;
        nCellsBx = 6;
        nCellsBy = 6;
        nCellsBz = 4;
        nAImpurity = 0;
        nAVacancy = 0;
        nBImpurity = 0;
        nBVacancy = 0;

        box = new Box(new BoundaryRectangularSlit(2, space), space);
        addBox(box);

        // beta-Sn box

        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the
        //lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815
        //angstroms) are taken from the ASM Handbook.
        aA = bA = 5.8314;
        cA = 3.1815;
        int basisA = 4;
        PrimitiveTetragonal primitiveA = new PrimitiveTetragonal(space, aA, cA);
        //Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
        //box.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
        //PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
        BravaisLatticeCrystal latticeA = new BravaisLatticeCrystal(primitiveA, new BasisBetaSnA5());

        aB = bB = 5.8314;
        cB = 3.1815;
        int basisB = 4;
        PrimitiveTetragonal primitiveB = new PrimitiveTetragonal(space, aB, cB);
        BravaisLatticeCrystal latticeB = new BravaisLatticeCrystal(primitiveB, new BasisBetaSnA5());

        box.getBoundary().setBoxSize(new Vector3D(aA * nCellsAx, aA * nCellsAy, (cA * nCellsAz) + (cB * nCellsBz)));


        //FCC Cu
        /**
         aA = bA = cA = 3.6148;SpeciesSpheresMono
         boxA.setDimensions(new Vector3D(aA*4, aA*4, aA*4));
         PrimitiveCubic primitiveA = new PrimitiveCubic(space, aA);
         LatticeCrystal latticeA = new LatticeCrystal(new Crystal(
         primitiveA, new BasisCubicFcc(primitiveA)));

         aB = bB = cB = 3.6148;
         boxB.setDimensions(new Vector3D(aB*4, aB*4, aB*4));
         PrimitiveCubic primitiveB = new PrimitiveCubic(space, aB);
         LatticeCrystal latticeB = new LatticeCrystal(new Crystal(
         primitiveB, new BasisCubicFcc(primitiveB)));

         */

        //FCC Ag
        /**
         aA = bA = cA = 4.0863;
         boxA.setDimensions(new Vector3D(aA*4, aA*4, aA*4));
         PrimitiveCubic primitiveA = new PrimitiveCubic(space, aA);
         LatticeCrystal latticeA = new LatticeCrystal(new Crystal(
         primitiveA, new BasisCubicFcc(primitiveA)));

         aB = bB = cB = 4.0863;
         boxB.setDimensions(new Vector3D(aB*4, aB*4, aB*4));
         PrimitiveCubic primitiveB = new PrimitiveCubic(space, aB);
         LatticeCrystal latticeB = new LatticeCrystal(new Crystal(
         primitiveB, new BasisCubicFcc(primitiveB)));
         */

        nA = (nCellsAx * nCellsAy * nCellsAz) * basisA;
        nAFixed = (nCellsAx * nCellsAy * 2) * basisA;
        nAMobile = nA - nAFixed - nAImpurity - nAVacancy;
        nB = (nCellsBx * nCellsBy * nCellsBz) * basisB;
        nBFixed = (nCellsBx * nCellsBy * 2) * basisB;
        nBMobile = nB - nBFixed - nBImpurity - nBVacancy;

        addSpecies(snFixedA);
        addSpecies(snA);
//        addSpecies(agA);
//        addSpecies(cuA);
        addSpecies(snFixedB);
        addSpecies(snB);
//        addSpecies(agB);
//        addSpecies(cuB);

        box.setNMolecules(snFixedA, nAFixed);
        box.setNMolecules(snA, nAMobile);
        box.setNMolecules(snFixedB, nBFixed);
        box.setNMolecules(snB, nBMobile);


        GrainBoundaryConfiguration config = new GrainBoundaryConfiguration(latticeA, latticeB, space);
        config.setDimensions(nCellsAx, nCellsAy, nCellsAz, nCellsBx, nCellsBy,
                nCellsBz, aA, bA, cA, aB, bB, cB);
        config.initializeCoordinates(box);

        potentialN = new PotentialMEAM(space);
        potentialN.setParameters(snFixedA.getLeafType(), ParameterSetMEAM.Sn);
        potentialN.setParameters(snA.getLeafType(), ParameterSetMEAM.Sn);
//		potentialN.setParameters(agA.getLeafType(), ParameterSetMEAM.Ag);
//		potentialN.setParameters(cuA.getLeafType(), ParameterSetMEAM.Cu);
        potentialN.setParameters(snFixedB.getLeafType(), ParameterSetMEAM.Sn);
        potentialN.setParameters(snB.getLeafType(), ParameterSetMEAM.Sn);
//		potentialN.setParameters(agB.getLeafType(), ParameterSetMEAM.Ag);
//		potentialN.setParameters(cuB.getLeafType(), ParameterSetMEAM.Cu);
//		potentialN.setParametersIMC(cuA.getLeafType(), ParameterSetMEAM.Cu3Sn);
//		potentialN.setParametersIMC(agA.getLeafType(), ParameterSetMEAM.Ag3Sn);
//		potentialN.setParametersIMC(cuB.getLeafType(), ParameterSetMEAM.Cu3Sn);
//		potentialN.setParametersIMC(agB.getLeafType(), ParameterSetMEAM.Ag3Sn);
//        this.potentialMaster.addPotential(potentialN, new Species[]{snFixedA, snA, agA, cuA, snFixedB, snB, agB, cuB});
        this.potentialMaster.addPotential(potentialN, new AtomType[]{snFixedA.getLeafType(), snA.getLeafType(), snFixedB.getLeafType(), snB.getLeafType()});
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
        MEAM_3DMDwithGB sim = new MEAM_3DMDwithGB();

        MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
        MeterKineticEnergy kineticMeter = new MeterKineticEnergy(sim.box);

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

        SimulationGraphic simgraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
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
        colorScheme.setColor(sim.snFixedA.getLeafType(), Color.blue);
        colorScheme.setColor(sim.snA.getLeafType(), Color.red);
//    	colorScheme.setColor(sim.agA.getMoleculeType(),java.awt.Color.gray);
//    	colorScheme.setColor(sim.cuA.getMoleculeType(),java.awt.Color.orange);
        colorScheme.setColor(sim.snFixedB.getLeafType(), Color.yellow);
        colorScheme.setColor(sim.snB.getLeafType(), Color.green);
//    	colorScheme.setColor(sim.agB.getMoleculeType(),java.awt.Color.gray);
//    	colorScheme.setColor(sim.cuB.getMoleculeType(),java.awt.Color.orange);

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
