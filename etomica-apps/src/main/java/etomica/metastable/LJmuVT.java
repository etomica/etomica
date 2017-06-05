/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.metastable;

import etomica.action.activity.ActivityIntegrate;
import etomica.integrator.IntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.statmech.LennardJones;
import etomica.units.Energy;
import etomica.units.Pixel;
import etomica.units.SimpleUnit;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.B3LJ;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class LJmuVT extends Simulation {
    
    public final PotentialMasterCell potentialMaster;
    public final SpeciesSpheresMono species;
    public final Box box;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;
    public final MCMoveInsertDelete mcMoveID;

    public LJmuVT(Space _space, int numAtoms, double temperature, double density, double mu, double rc, int[] seeds) {
        super(_space);
        if (seeds != null) {
            setRandom(new RandomMersenneTwister(seeds));
        }
        potentialMaster = new PotentialMasterCell(this, rc, space);
        potentialMaster.setCellRange(2);

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
        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});

        //construct box
	    box = new Box(space);
        addBox(box);
        integrator.setBox(box);

        mcMoveID = new MCMoveInsertDelete(potentialMaster, random, space);
        mcMoveID.setMu(mu);
        mcMoveID.setSpecies(species);
        integrator.getMoveManager().addMCMove(mcMoveID);
        
        double L = Math.pow(numAtoms/density, 1.0/3.0);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{L,L,L}));
        box.setNMolecules(species, numAtoms);
        
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);

        potentialMaster.getNbrCellManager(box).assignCellAll();
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }
    
    public static double[] calcBMu(double density, double B2, double B3) {
        double muCalc = Math.log(density);
        double dmuCalc = 1/density;
        muCalc += 2 * B2 * density;
        dmuCalc += 2 * B2;
        muCalc += 1.5 * B3 * density * density;
        dmuCalc += 3 * B3 * density;
        return new double[]{muCalc,dmuCalc};
    }

    public static void main(String[] args) throws IOException {
        LJMCParams params = new LJMCParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.mu = -2.8;
        }
        Space space = Space3D.getInstance();
        final int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double mu = params.mu;
        double mu0 = params.mu0;
        long numSteps = params.numSteps;
        int numRuns = params.numRuns;
        double rc = params.rc;

        double B2 = LennardJones.B2(temperature);
        double B3 = B3LJ.value(temperature)[0];

        double bmu = mu/temperature;
        double density1 = Math.exp(bmu);
        double density2 = 2*density1;
        double[] bdmu2 = calcBMu(density2, B2, B3);
        double bmu2 = bdmu2[0];
        boolean fudged = false;
        while (bdmu2[1] < 0) {
            fudged = true;
            bmu -= 0.01;
            density1 = Math.exp(bmu);
            density2 = 2*density1;
            bdmu2 = calcBMu(density2, B2, B3);
            bmu2 = bdmu2[0];
        }
        while (bmu2 < bmu) {
            density2 += 0.01*density1;
            bdmu2 = calcBMu(density2, B2, B3);
            bmu2 = bdmu2[0];
            while (bdmu2[1] < 0) {
                fudged = true;
                bmu -= 0.01;
                density1 = Math.exp(bmu);
                density2 = 2*density1;
                bdmu2 = calcBMu(density2, B2, B3);
                bmu2 = bdmu2[0];
            }

        }
        if (fudged) {
            System.out.println("decreased initial chemical potential estimate to "+bmu*temperature);
        }
        double density = 0;
        int iter = 0;
        while (true) {
            iter++;
            double density3 = 0.5*(density1 + density2);
            double bmu3 = calcBMu(density3, B2, B3)[0];
            if (Math.abs(bmu3-bmu) < 0.01) {
                density = density3;
                break;
            }
            if (bmu3 > bmu) {
                density2 = density3;
            }
            else {
                density1 = density3;
            }
            if (iter>100) {
                bmu -= 0.01;
                throw new RuntimeException("unable to estimate density");
            }
        }
        double vol = numAtoms/density;
        if (Math.pow(vol,1.0/3.0) < 2*rc) {
            System.out.println("rc: "+rc);
            System.out.println("density: "+density+"   volume: "+vol);
            throw new RuntimeException("System size seems too small for cutoff");
        }
        int numAtoms0 = numAtoms/10;
        double density0 = numAtoms0/vol;


        final FileWriter fw = params.out == null ? null : new FileWriter(params.out);
        
        System.out.println("N: "+numAtoms);
        System.out.println("mu: "+mu);
        System.out.println("numSteps: "+numSteps);
        System.out.println("rc: "+rc);
        System.out.println("T: "+temperature);
        System.out.println("V0: "+vol);

        for (int i=0; i<numRuns; i++){

            int[] seeds = null; //new int[]{-475498437, 1044754174, -1894223345, 1180064439};
            final LJmuVT sim = new LJmuVT(space, numAtoms0, temperature, density0, mu0, rc, seeds);
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
            
            final MeterNMolecules meterN = new MeterNMolecules();
            meterN.setBox(sim.box);

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
                avgDensity.addDataSink(historyDensity, new StatType[]{AccumulatorAverage.MOST_RECENT});
                DisplayPlot densityPlot = new DisplayPlot();
                densityPlot.setLabel("density");
                historyDensity.addDataSink(densityPlot.getDataSet().makeDataSink());
                ljmcGraphic.add(densityPlot);
                densityPlot.setDoLegend(false);

                AccumulatorHistory historyN = new AccumulatorHistory(new HistoryCollapsingAverage());
                DataPumpListener pumpN = new DataPumpListener(meterN, historyN, numAtoms);
                sim.integrator.getEventManager().addListener(pumpN);
                dataStreamPumps.add(pumpN);
                DisplayPlot nPlot = new DisplayPlot();
                nPlot.setLabel("N");
                historyN.addDataSink(nPlot.getDataSet().makeDataSink());
                ljmcGraphic.add(nPlot);
                nPlot.setDoLegend(false);

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
                
                final DeviceSlider muSlider = new DeviceSlider(sim.getController(), sim.mcMoveID, "mu");
                muSlider.setPrecision(2);
                muSlider.setNMajor(5);
                muSlider.setMinimum(-5);
                muSlider.setMaximum(0);
                muSlider.setValue(mu0);
                muSlider.setShowValues(true);
                muSlider.setEditValues(true);
                muSlider.setShowBorder(true);
                muSlider.setLabel("mu");
                ljmcGraphic.add(muSlider);

                return;
            }
            long tstart = System.currentTimeMillis();
            sim.activityIntegrate.setMaxSteps(numSteps/10);
            sim.getController().actionPerformed();
            sim.integrator.resetStepCount();
            sim.getController().reset();
            if (numRuns==1) System.out.println("Equilibration finished");

            sim.mcMoveID.setMu(mu);
            sim.activityIntegrate.setMaxSteps(numSteps);
            sim.integrator.getEventManager().addListener(new IIntegratorListener() {

                int count = 0;
                int interval = numAtoms;
                
                public void integratorStepStarted(IntegratorEvent e) {
                }
                
                public void integratorStepFinished(IntegratorEvent e) {
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
                
                public void integratorInitialized(IntegratorEvent e) {
                }
            });
            
            sim.getController().actionPerformed();
            long tstop = System.currentTimeMillis();
            
            if (numRuns <= 10 || (numRuns*10/(i+1))*(i+1) == numRuns*10) {
                System.out.println("Run "+(i+1)+" finished");
            }
            System.out.println("time: "+(tstop-tstart)*0.001);
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
        public int numAtoms = 400;
        public double temperature = 0.7;
        public double mu0 = -4.3;
        public double mu = -3;
        public long numSteps = 10000000;
        public int numRuns = 10;
        public double rc = 4;
        public String out = "rhoPU.dat";
    }
}
