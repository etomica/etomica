/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.metastable;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtomType;
import etomica.box.Box;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataDump;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Energy;
import etomica.units.Pixel;
import etomica.units.SimpleUnit;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.RandomMersenneTwister;

public class LJNPT extends Simulation {
    
    public final PotentialMasterCell potentialMaster;
    public final SpeciesSpheresMono species;
    public final Box box;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;
    public final MCMoveAtom mcMoveAtom;
    public final MCMoveVolume mcMoveVolume;

    public LJNPT(Space _space, int numAtoms, double temperature, double density, double pressure, double rc, int[] seeds) {
        super(_space);
        if (seeds != null) {
            setRandom(new RandomMersenneTwister(seeds));
        }
        potentialMaster = new PotentialMasterCell(this, rc, space);
        potentialMaster.setCellRange(2);
        potentialMaster.lrcMaster().setEnabled(false);
        
        //controller and integrator
	    integrator = new IntegratorMC(potentialMaster, random, temperature);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

	    //species and potentials
	    species = new SpeciesSpheresMono(this, space);//index 1
	    species.setIsDynamic(true);
        addSpecies(species);
        
        //instantiate several potentials for selection in combo-box
	    P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential, rc);
	    potentialMaster.addPotential(p2Truncated, new IAtomType[]{species.getLeafType(), species.getLeafType()});
	    
        //construct box
	    box = new Box(space);
        addBox(box);
        integrator.setBox(box);

        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        
        mcMoveVolume = new MCMoveVolume(potentialMaster, random, space, pressure);
        mcMoveVolume.setTemperature(temperature);
        integrator.getMoveManager().addMCMove(mcMoveVolume);
        integrator.getMoveManager().setFrequency(mcMoveVolume, 5);
        
        double L = Math.pow(numAtoms/density, 1.0/3.0);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{L,L,L}));
        box.setNMolecules(species, numAtoms);
        
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);

        potentialMaster.getNbrCellManager(box).assignCellAll();
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }
    
    public static void main(String[] args) throws IOException {
        LJMCParams params = new LJMCParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
        }
        Space space = Space3D.getInstance();
        final int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double density = params.density;
        double pressure0 = params.pressure0;
        double pressure = params.pressure;
        long numSteps = params.numSteps;
        int numRuns = params.numRuns;
        double rc = params.rc;

        final FileWriter fw = params.out == null ? null : new FileWriter(params.out);
        
        System.out.println("N: "+numAtoms);
        System.out.println("numSteps: "+numSteps);
        System.out.println("rc: "+rc);
        System.out.println("T: "+temperature);
        System.out.println("P: "+pressure);

        for (int i=0; i<numRuns; i++){

            int[] seeds = null; //new int[]{-475498437, 1044754174, -1894223345, 1180064439};
            final LJNPT sim = new LJNPT(space, numAtoms, temperature, density, pressure0, rc, seeds);
            if (seeds == null) {
                System.out.println("random seeds: "+Arrays.toString(sim.getRandomSeeds()));
            }
            else {
                System.out.println("set random seeds: "+Arrays.toString(seeds));
            }
            
            final MeterPressure meterP = new MeterPressure(space);
            meterP.setBox(sim.box);
            meterP.setPotentialMaster(sim.potentialMaster);
            meterP.setTemperature(temperature);
    
            final MeterDensity meterDensity = new MeterDensity(sim.getSpace());
            meterDensity.setBox(sim.box);
    
            final MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator();
            meterPE.setIntegrator(sim.integrator);
            
            if (false) {
                SimulationGraphic ljmcGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
                ljmcGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(3));
                
                SimulationGraphic.makeAndDisplayFrame(ljmcGraphic.getPanel(), "LJ Spinodal");
                
                List<DataPump> dataStreamPumps = ljmcGraphic.getController().getDataStreamPumps();
                
                AccumulatorAverageCollapsing avgDensity = new AccumulatorAverageCollapsing();
                avgDensity.setPushInterval(1);
                DataPumpListener pumpDensity = new DataPumpListener(meterDensity, avgDensity, numAtoms);
                sim.integrator.getEventManager().addListener(pumpDensity);
                dataStreamPumps.add(pumpDensity);
                DisplayTextBoxesCAE displayDensity = new DisplayTextBoxesCAE();
                displayDensity.setAccumulator(avgDensity);
                ljmcGraphic.add(displayDensity);
                
                
                AccumulatorHistory historyDensity = new AccumulatorHistory(new HistoryCollapsingAverage());
                avgDensity.addDataSink(historyDensity, new StatType[]{avgDensity.MOST_RECENT});
                DisplayPlot densityPlot = new DisplayPlot();
                densityPlot.setLabel("density");
                historyDensity.addDataSink(densityPlot.getDataSet().makeDataSink());
                ljmcGraphic.add(densityPlot);
                densityPlot.setDoLegend(false);
                
                DataDump densityDump = new DataDump();
                historyDensity.addDataSink(densityDump);
                
                AccumulatorHistory historyEnergy = new AccumulatorHistory(new HistoryCollapsingAverage());
                DataPumpListener pumpEnergy = new DataPumpListener(meterPE, historyEnergy, numAtoms);
                sim.integrator.getEventManager().addListener(pumpEnergy);
                dataStreamPumps.add(pumpEnergy);
                DisplayPlot energyPlot = new DisplayPlot();
                energyPlot.setUnit(new SimpleUnit(Energy.DIMENSION, numAtoms, "energy", "U", false));
                energyPlot.setLabel("PE");
                historyEnergy.addDataSink(energyPlot.getDataSet().makeDataSink());
                ljmcGraphic.add(energyPlot);
                energyPlot.setDoLegend(false);
                
                DataProcessorXY dpXY = new DataProcessorXY(densityDump);
                historyEnergy.addDataSink(dpXY);
                DisplayPlot densityEnergyPlot = new DisplayPlot();
                dpXY.setDataSink(densityEnergyPlot.getDataSet().makeDataSink());
                densityEnergyPlot.setLabel("Density/Energy");
                densityEnergyPlot.getPlot().setYLabel("Potential Energy");
                densityEnergyPlot.getPlot().setXLabel("Density");
                densityEnergyPlot.setUnit(new SimpleUnit(Energy.DIMENSION, numAtoms, "energy", "U", false));
                ljmcGraphic.add(densityEnergyPlot);
                
                AccumulatorHistory historyPressure = new AccumulatorHistory(new HistoryCollapsingAverage());
                DataPumpListener pumpPressure = new DataPumpListener(meterP, historyPressure, numAtoms);
                sim.integrator.getEventManager().addListener(pumpPressure);
                dataStreamPumps.add(pumpPressure);
                DisplayPlot pressurePlot = new DisplayPlot();
                pressurePlot.setLabel("Pressure");
                historyPressure.addDataSink(pressurePlot.getDataSet().makeDataSink());
                ljmcGraphic.add(pressurePlot);
                pressurePlot.setDoLegend(false);
                
                final DeviceSlider pSlider = new DeviceSlider(sim.getController(), sim.mcMoveVolume, "pressure");
                pSlider.setPrecision(4);
                pSlider.setNMajor(4);
                pSlider.setMaximum(0.1);
                pSlider.setValue(pressure0);
                pSlider.setShowValues(true);
                pSlider.setEditValues(true);
                pSlider.setShowBorder(true);
                pSlider.setLabel("Pressure");
                ljmcGraphic.add(pSlider);
                
                return;
            }
            
            sim.activityIntegrate.setMaxSteps(numSteps/10);
            sim.getController().actionPerformed();
            sim.integrator.resetStepCount();
            sim.getController().reset();
            if (numRuns==1) System.out.println("Equilibration finished");

            sim.mcMoveVolume.setPressure(pressure);
            sim.activityIntegrate.setMaxSteps(numSteps);
            sim.integrator.getEventManager().addListener(new IIntegratorListener() {

                int count = 0;
                int interval = numAtoms;
                
                public void integratorStepStarted(IIntegratorEvent e) {
                }
                
                public void integratorStepFinished(IIntegratorEvent e) {
                    interval--;
                    if (interval > 0) return;
                    interval = numAtoms;
                    double U = meterPE.getDataAsScalar()/numAtoms;
                    double P = meterP.getDataAsScalar();
                    double rho = meterDensity.getDataAsScalar();
                    try {
                        if (fw!=null) fw.write(count+" "+rho+" "+P+" "+U+"\n");
                        else System.out.println(count+" "+rho+" "+P+" "+U);
                    }
                    catch (IOException ex) {
                        throw new RuntimeException(ex);
                    }
                    count++;
                    if (U<-1) sim.activityIntegrate.setMaxSteps(0);
                }
                
                public void integratorInitialized(IIntegratorEvent e) {
                }
            });
            
            sim.getController().actionPerformed();
            
            if (numRuns <= 10 || (numRuns*10/(i+1))*(i+1) == numRuns*10) {
                System.out.println("Run "+(i+1)+" finished");
            }
            if (fw != null && numRuns > 1) {
                fw.write("\n");
            }
            if (seeds != null) break;
        }

        if (fw != null) {
            fw.close();
        }
    }
    
    public static class LJMCParams extends ParameterBase {
        public int numAtoms = 2000;
        public double temperature = 0.7;
        public double density = 0.002;
        public double pressure0 = 0.001;
        public double pressure = 0.01;
        public long numSteps = 10000000;
        public int numRuns = 10;
        public double rc = 4;
        public String out = "rhoPU.dat";
    }
}
