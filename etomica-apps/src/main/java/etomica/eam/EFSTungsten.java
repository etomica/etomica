/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.eam;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.Tungsten;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.history.History;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.graphics.*;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.crystal.BasisCubicBcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.CriterionSimple;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.MeterDADB;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.*;
import etomica.units.dimensions.Energy;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Molecular-Dynamics Simulation Using the extended Finnis-Sinclair Method 
 * (EFS) Potential.  
 * 
 * The EFS potential is intended for use with metallic and covalently-bonded
 * solid systems.
 * 
 * The EFS potential for an atom is built using terms describing parts of the
 * relationships between the atom and each of its neighbors, the number of which 
 * is determined by a cutoff variable.  Each type of pair-
 * wise term is summed over all the neighbors, and then used in expressions 
 * describing the embedding energy and the repulsive energy of the atom.   
 * Effectively, the EFS potential is a many-body potential.  
 * 
 * This class was adapted from LjMd3D.java by K.R. Schadel and A. Schultz in July 
 * 2005.  Intitially, it employed a version of the embedded-atom method potential, 
 * and was later adapted in February 2006 to use the modified embedded-atom method
 * potential. Finally it was modified to use the extended Finnis-Sinclair Method in 
 * October 2013.
 * 
 * @author Joe Kromph
 */
 
public class EFSTungsten extends Simulation {
    
    private static final String APP_NAME = "MEAM Md3D";
    public final PotentialMasterList potentialMaster;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono w;
    public Box box;
    public PotentialEFS potentialN;
    public Controller controller;
    public DisplayBox display;
    public MeterEnergy energy;
    public ActivityIntegrate activityIntegrate;
    public IDataInfo info2;
    public CoordinateDefinition coordinateDefinition;

    public EFSTungsten(int numatoms, double density, double temperature) {
        super(Space3D.getInstance());

        w = new SpeciesSpheresMono(space, Tungsten.INSTANCE);
        w.setIsDynamic(true);
        addSpecies(w);

        potentialMaster = new PotentialMasterList(this, space);
        box = this.makeBox();
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setTimeStep(0.001);
        integrator.setTemperature(Kelvin.UNIT.toSim(temperature));
        integrator.setThermostatInterval(100);
        integrator.setIsothermal(true);
        integrator.setThermostatNoDrift(true);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        box.setNMolecules(w, numatoms);

        //BCC W

        BoxInflate inflate = new BoxInflate(box, space);
        inflate.setTargetDensity(density); //(atoms/A^3)
        inflate.actionPerformed();
        int nCells = (int) Math.round(Math.pow(numatoms / 2, 1.0 / 3.0));
        coordinateDefinition = new CoordinateDefinitionLeaf(box, new PrimitiveCubic(space, box.getBoundary().getBoxSize().getX(0) / nCells), new BasisCubicBcc(), space);
        coordinateDefinition.initializeCoordinates(new int[]{nCells, nCells, nCells});


        double c, c0, c1, c2, c3, c4, A, d, B;

        c = 3.25;
        c0 = ElectronVolt.UNIT.toSim(48.52796);
        c1 = ElectronVolt.UNIT.toSim(-33.79621);
        c2 = ElectronVolt.UNIT.toSim(5.854334);
        c3 = ElectronVolt.UNIT.toSim(-0.0098221);
        c4 = ElectronVolt.UNIT.toSim(0.033338);
        A = ElectronVolt.UNIT.toSim(1.885948);
        d = 4.41;
        B = 0;

        potentialN = new PotentialEFS(space, A, B, c, d, c0, c1, c2, c3, c4);

        this.potentialMaster.addPotential(potentialN, new AtomType[]{w.getLeafType()});
        potentialMaster.setRange(potentialN.getRange() * 1.3);
        potentialMaster.setCriterion(potentialN, new CriterionSimple(this, space, potentialN.getRange(), potentialN.getRange() * 1.3));
//        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        potentialMaster.getNeighborManager(box).reset();
    }

    public static void main(String[] args) {
        Params parameters=new Params();
        ParseArgs.doParseArgs(parameters, args);
        int numatoms=parameters.numatoms;
        double density=parameters.density;
        //System.out.println(Mole.UNIT.fromSim(16.870e24));
        double temperature=parameters.temperature;
        long numsteps=parameters.numsteps;
        boolean doHistory = parameters.doHistory;
        
    	EFSTungsten sim = new EFSTungsten(numatoms, density, temperature);
    	
    	System.out.println("EFS Tungsten at rho="+density+" and T="+temperature+"K");
    	System.out.println(numatoms+" atoms");
    	System.out.println(numsteps+" steps");
    	
    	MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
        MeterKineticEnergy kineticMeter = new MeterKineticEnergy(sim.box);
        MeterEnergy totalEnergyMeter = new MeterEnergy(sim.potentialMaster,sim.box);
        MeterDADB DADBMeter = new MeterDADB(sim.getSpace(),energyMeter,sim.potentialMaster,sim.coordinateDefinition,Kelvin.UNIT.toSim(temperature));
        MeterPressure pressureMeter = new MeterPressure(sim.getSpace());
    	System.out.println("lattice energy "+ElectronVolt.UNIT.fromSim(energyMeter.getDataAsScalar()/numatoms));
    	System.out.println("eV is "+ElectronVolt.UNIT.fromSim(1));
    	System.out.println("K is "+Kelvin.UNIT.fromSim(1));
        DADBMeter.setLatticeEnergy(energyMeter.getDataAsScalar());
        MeterDADB.justU = true;
        pressureMeter.setBox(sim.box);
        pressureMeter.setIntegrator(sim.integrator);
        DADBMeter.getData();
                
        AccumulatorHistory energyAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorHistory kineticAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorHistory totalEnergyAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorHistory DADBAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorHistory pressureAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        
        
    	if(false){ //graph?
            AccumulatorAverageCollapsing accumulatorAveragePE = new AccumulatorAverageCollapsing();
            AccumulatorAverageCollapsing accumulatorAverageKE = new AccumulatorAverageCollapsing();
            AccumulatorAverageCollapsing accumulatorAverageDADB = new AccumulatorAverageCollapsing();
            AccumulatorAverageCollapsing accumulatorAveragePressure = new AccumulatorAverageCollapsing();
            
            DataPumpListener energyPump = new DataPumpListener(energyMeter,accumulatorAveragePE,10);    
            DataPumpListener kineticPump = new DataPumpListener(kineticMeter, accumulatorAverageKE);
            DataPumpListener totalEnergyPump = new DataPumpListener(totalEnergyMeter, totalEnergyAccumulator,10);
            DataPumpListener DADBPump = new DataPumpListener(DADBMeter, DADBAccumulator,10);
            DataPumpListener pressurePump = new DataPumpListener(pressureMeter, pressureAccumulator,10);
            
            sim.integrator.getEventManager().addListener((energyPump));
            sim.integrator.getEventManager().addListener((kineticPump));
            sim.integrator.getEventManager().addListener((totalEnergyPump));
            sim.integrator.getEventManager().addListener((DADBPump));
            sim.integrator.getEventManager().addListener((pressurePump));

            accumulatorAveragePE.addDataSink(energyAccumulator, new StatType[]{AccumulatorAverage.MOST_RECENT});
            accumulatorAverageKE.addDataSink(kineticAccumulator, new StatType[]{AccumulatorAverage.MOST_RECENT});
            accumulatorAverageDADB.addDataSink(DADBAccumulator, new StatType[]{AccumulatorAverage.MOST_RECENT});
            accumulatorAveragePressure.addDataSink(pressureAccumulator, new StatType[]{AccumulatorAverage.MOST_RECENT});

            DisplayPlot plotPE = new DisplayPlot();
            plotPE.setLabel("PE Plot");
            plotPE.setUnit(new SimpleUnit(Energy.DIMENSION, numatoms, "", "", false));
            DisplayPlot plotKE = new DisplayPlot();
            plotKE.setLabel("KE Plot");
            DisplayPlot plottotalEnergy = new DisplayPlot();
            plottotalEnergy.setLabel("Total Energy Plot");
            DisplayPlot plotPressure = new DisplayPlot();
            plotPressure.setLabel("Pressure");
        	
            energyAccumulator.setDataSink(plotPE.getDataSet().makeDataSink());
            kineticAccumulator.setDataSink(plotKE.getDataSet().makeDataSink());
            totalEnergyAccumulator.setDataSink(plottotalEnergy.getDataSet().makeDataSink());
            DADBAccumulator.setDataSink(plotPE.getDataSet().makeDataSink());
            pressureAccumulator.setDataSink(plotPressure.getDataSet().makeDataSink());
            
            //energyAccumulator.setBlockSize(50);
        	accumulatorAveragePE.setPushInterval(1);
        	accumulatorAverageKE.setPushInterval(1);
        	accumulatorAverageDADB.setPushInterval(1);
        	accumulatorAveragePressure.setPushInterval(1);
        	
        	//Heat Capacity (PE)
        	DataProcessorCvMD dataProcessorPE = new DataProcessorCvMD();
        	dataProcessorPE.setIntegrator(sim.integrator);
        	
        	//Heat Capacity (KE)
        	DataProcessorCvMD dataProcessorKE = new DataProcessorCvMD();
        	dataProcessorKE.setIntegrator(sim.integrator);

            accumulatorAveragePE.addDataSink(dataProcessorPE, new StatType[]{AccumulatorAverage.STANDARD_DEVIATION});
            accumulatorAverageKE.addDataSink(dataProcessorKE, new StatType[]{AccumulatorAverage.STANDARD_DEVIATION});
           
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
            ArrayList<DataPump> dataStreamPumps = simGraphic.getController().getDataStreamPumps();
            dataStreamPumps.add(energyPump);
            dataStreamPumps.add(kineticPump);
            dataStreamPumps.add(DADBPump);
            dataStreamPumps.add(pressurePump);
            
        	DisplayTextBox cvBoxPE = new DisplayTextBox();
        	dataProcessorPE.setDataSink(cvBoxPE);
        	cvBoxPE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double []{1,-1,-1}));
        	cvBoxPE.setLabel("PE Cv contrib.");
        	DisplayTextBox cvBoxKE = new DisplayTextBox();
        	dataProcessorKE.setDataSink(cvBoxKE);
        	cvBoxKE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double []{1,-1,-1}));
        	cvBoxKE.setLabel("KE Cv contrib.");
        	plotPressure.setUnit(Bar.UNIT);
    
        	simGraphic.add(plotPE);
        	simGraphic.add(plotKE);
            simGraphic.add(plottotalEnergy);
            simGraphic.add(plotPressure);
                    	
        	simGraphic.getPanel().controlPanel.add(cvBoxKE.graphic(), SimulationPanel.getVertGBC());
        	simGraphic.getPanel().controlPanel.add(cvBoxPE.graphic(), SimulationPanel.getVertGBC());
    
        	simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
    
        	simGraphic.makeAndDisplayFrame(APP_NAME);
        	return;
    	}
    	
    	sim.activityIntegrate.setMaxSteps(numsteps/10);
    	sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numsteps);
        System.out.println("equilibration finished");

        AccumulatorAverageFixed accumulatorAveragePE = null;
        AccumulatorAverageFixed accumulatorAverageDADB = null;
        AccumulatorHistory accumulatorHistoryPE = null;
        AccumulatorHistory accumulatorHistoryDADB = null;

        accumulatorAveragePE = new AccumulatorAverageFixed((numsteps/10)/100);
        accumulatorAverageDADB = new AccumulatorAverageFixed((numsteps/10)/100);

        DataPumpListener energyPump = new DataPumpListener(energyMeter,accumulatorAveragePE,10);
        DataPumpListener dadbPump = new DataPumpListener(DADBMeter, accumulatorAverageDADB, 10);

        sim.integrator.getEventManager().addListener(energyPump);
        sim.integrator.getEventManager().addListener(dadbPump);
        if (doHistory) {
            accumulatorHistoryPE = new AccumulatorHistory(new HistoryCollapsingAverage(200));
            accumulatorHistoryDADB = new AccumulatorHistory(new HistoryCollapsingAverage(200));
            DataSourceCountTime timer = new DataSourceCountTime(sim.integrator);
            accumulatorHistoryPE.setTimeDataSource(timer);
            accumulatorHistoryDADB.setTimeDataSource(timer);

            energyPump = new DataPumpListener(energyMeter,accumulatorHistoryPE,10);
	        dadbPump = new DataPumpListener(DADBMeter, accumulatorHistoryDADB, 10);

	        sim.integrator.getEventManager().addListener(energyPump);
	        sim.integrator.getEventManager().addListener(dadbPump);        	
        }
        
        
        long t1 = System.currentTimeMillis();
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

        if (doHistory) {
        	History histPE = accumulatorHistoryPE.getHistory();
        	History histDADB = accumulatorHistoryDADB.getHistory();
        	double[] x = histPE.getXValues();
        	double[] hpe = histPE.getHistory();
        	double[] hdadb = histDADB.getHistory();
        	try {
        	    double fac = ElectronVolt.UNIT.fromSim(1.0/(numatoms));
        		FileWriter hWriter = new FileWriter("pe.dat");
        		for (int i=0; i<x.length; i++) {
        		    if (!Double.isNaN(hpe[i])) {
        		        hWriter.write(x[i]+" "+(hpe[i]+1.5*temperature)*fac+"\n");
        		    }
        		}
        		hWriter.close();
                hWriter = new FileWriter("dadb.dat");
                for (int i=0; i<x.length; i++) {
                    if (!Double.isNaN(hdadb[i])) {
                        hWriter.write(x[i]+" "+hdadb[i]*fac+"\n");
                    }
                }
                hWriter.close();
        	}
        	catch (IOException e) {
        		throw new RuntimeException(e);
        	}
        }

        double PE = accumulatorAveragePE.getData(AccumulatorAverage.AVERAGE).getValue(0)
                    /sim.box.getLeafList().size();
        double PEerror = accumulatorAveragePE.getData(AccumulatorAverage.ERROR).getValue(0)
                /sim.box.getLeafList().size();
        double PEcor = accumulatorAveragePE.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0)
                /sim.box.getLeafList().size();
        double PV = Pascal.UNIT.toSim(0)*1e9/(Mole.UNIT.toSim(1/16.870e24));

        double dadb = accumulatorAverageDADB.getData(AccumulatorAverage.AVERAGE).getValue(0)
                /sim.box.getLeafList().size();
        double dadbError = accumulatorAverageDADB.getData(AccumulatorAverage.ERROR).getValue(0)
                /sim.box.getLeafList().size();
        double dadbCor = accumulatorAverageDADB.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0)
                /sim.box.getLeafList().size();
	
        System.out.println("PE(eV) "+ElectronVolt.UNIT.fromSim(PE)+" error: "+ElectronVolt.UNIT.fromSim(PEerror)+ " corrolation: "+PEcor);
	//        System.out.println("PV(ev)= "+ElectronVolt.UNIT.fromSim(PV));
        System.out.println("dadb(eV) "+ElectronVolt.UNIT.fromSim(dadb)+" error: "+ElectronVolt.UNIT.fromSim(dadbError)+ " corrolation: "+dadbCor);

        System.out.println("time: "+(t2-t1)*0.001);
    }

    public static class Params extends ParameterBase {
        public int numatoms = 2 * 4 * 4 * 4;
        public double density = 1.0 / 16.870; //2/3.151/3.151/3.151;//2/3.15961/3.15961/3.15961;//Mole.UNIT.toSim(1/16.870e24);
        public double temperature = 300;
        public long numsteps = 10000;
        public boolean doHistory = true;
    }
}
