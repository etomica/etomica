/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.histogram.Histogram;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterRDFNeighbors;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class MappedVirialLJ extends Simulation {

    public SpeciesGeneral species;
    public Box box;
    public IntegratorMC integrator;
    public MCMoveAtom move;
    public PotentialMasterCell potentialMaster;

    public P2SoftSphericalTruncated p2Truncated;

    public MappedVirialLJ(Space _space, int numAtoms, double temperature, double density, double rc) {
        super(_space);

        //species
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);
        box = this.makeBox();

        potentialMaster = new PotentialMasterCell(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        potentialMaster.doAllTruncationCorrection = false;

        //controller and integrator
        integrator = new IntegratorMC(potentialMaster, random, temperature, box);
        this.getController().addActivity(new ActivityIntegrate(integrator));
        move = new MCMoveAtom(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(move);

        //potentials
        P2LennardJones potential = new P2LennardJones();
        p2Truncated = new P2SoftSphericalTruncated(potential, rc);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), p2Truncated);

        //construct box
        box.setNMolecules(species, numAtoms);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
    }
    
    public static void main(String[] args) throws IOException {
        
        LJMDParams params = new LJMDParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.temperature = 2;
            params.density = 0.10;
            params.numSteps = 10000000;
            params.rc = 4;
            params.numAtoms = 200;
            params.functionsFile = "0.10";
        }

        int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double density = params.density;
        long numSteps = params.numSteps;
        double rc = params.rc;
        String functionsFile = params.functionsFile;
        boolean computeP = params.computeP;
        boolean computePMA = params.computePMA;
        boolean graphics = params.graphics;
        int nBlocks = params.nBlocks;

        System.out.println("Virial mapped average");
        System.out.println(numAtoms+" atoms, "+numSteps+" steps");
        System.out.println("density: "+density);
        System.out.println("temperature: "+temperature);
        System.out.println("cutoff: "+rc);
        System.out.println(nBlocks+" blocks");

        Space space = Space.getInstance(3);

        MappedVirialLJ sim = new MappedVirialLJ(space, numAtoms, temperature, density, rc);
        
        if (graphics) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.makeAndDisplayFrame();

            ArrayList<DataPump> dataStreamPumps = simGraphic.getController().getDataStreamPumps();
            
            MeterRDFNeighbors meterRDF = new MeterRDFNeighbors(sim.box, sim.potentialMaster.getCellManager());
            int nbins = (int)Math.round(rc/0.01);
            meterRDF.getXDataSource().setNValues(nbins);
            meterRDF.getXDataSource().setXMax(rc);

            DisplayPlot rdfPlot = new DisplayPlot();
            DataPump rdfPump = new DataPump(meterRDF,rdfPlot.getDataSet().makeDataSink());
            IntegratorListenerAction rdfPumpListener = new IntegratorListenerAction(rdfPump);
            sim.integrator.getEventManager().addListener(rdfPumpListener);
            rdfPumpListener.setInterval(numAtoms);
            dataStreamPumps.add(rdfPump);
            
            rdfPlot.setDoLegend(false);
            rdfPlot.getPlot().setTitle("Radial Distribution Function");
            rdfPlot.setLabel("RDF");
            simGraphic.add(rdfPlot);

            return;
        }
        
        long t1 = System.currentTimeMillis();

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps / 10));

        sim.integrator.getMoveManager().setEquilibrating(false);

        int nBins = 1000000;
        long numSamples = numSteps/numAtoms;
        long samplesPerBlock = numSamples/nBlocks;
        if (samplesPerBlock == 0) samplesPerBlock = 1;

        final MeterMappedVirial meterMappedVirial = new MeterMappedVirial(sim.integrator.getPotentialCompute(), sim.box, nBins);
        meterMappedVirial.getPotentialCallback().setTemperature(sim.integrator.getTemperature(), sim.p2Truncated);
        final AccumulatorAverageFixed accMappedVirial = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpMappedVirial = new DataPumpListener(meterMappedVirial, accMappedVirial, numAtoms);
        if (computePMA) sim.integrator.getEventManager().addListener(pumpMappedVirial);
        System.out.println("x0: "+meterMappedVirial.getPotentialCallback().getX0());
        double qu = meterMappedVirial.getPotentialCallback().getQU();
        System.out.println("qu: "+qu);
        double q = meterMappedVirial.getPotentialCallback().getQ();
        System.out.println("q: "+q);

        final MeterPressure meterP = new MeterPressure(sim.box, sim.potentialMaster);
        meterP.setTemperature(temperature);
        final AccumulatorAverageFixed accP = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpP = new DataPumpListener(meterP, accP, numAtoms);
        if (computeP) sim.integrator.getEventManager().addListener(pumpP);

        MeterMeanForce meterF = null;
        MeterRDFNeighbors meterRDF = null;
        if (functionsFile != null) {
            int nbins = (int)Math.round(rc/0.01);
            meterF = new MeterMeanForce(sim.integrator.getPotentialCompute(), sim.p2Truncated, sim.box, nbins);
            if (computeP && computePMA) sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterF, numAtoms));

            meterRDF = new MeterRDFNeighbors(sim.box, sim.potentialMaster.getCellManager());
            meterRDF.getXDataSource().setNValues(nbins);
            meterRDF.getXDataSource().setXMax(rc);
            sim.integrator.getEventManager().addListener(new DataPumpListener(meterRDF, null, numAtoms));
        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        IData mappedAvg = accMappedVirial.getData(accMappedVirial.AVERAGE);
        IData mappedErr = accMappedVirial.getData(accMappedVirial.ERROR);
        IData mappedCor = accMappedVirial.getData(accMappedVirial.BLOCK_CORRELATION);
        double avg = mappedAvg.getValue(0);
        double err = mappedErr.getValue(0);
        double cor = mappedCor.getValue(0);
        if (computePMA) System.out.print(String.format("avg: %13.6e   err: %11.4e   cor: % 4.2f\n", avg, err, cor));

        double pAvg = accP.getData(accP.AVERAGE).getValue(0);
        double pErr = accP.getData(accP.ERROR).getValue(0);
        double pCor = accP.getData(accP.BLOCK_CORRELATION).getValue(0);
        if (computeP) System.out.print(String.format("Pressure     avg: %13.6e   err: %11.4e   cor: % 4.2f\n", pAvg, pErr, pCor));

        if (functionsFile != null) {
            FileWriter fw = new FileWriter(functionsFile+"_gr.dat");
            IData rdata = meterRDF.getIndependentData(0);
            IData gdata = meterRDF.getData();
            for (int i=0; i<rdata.getLength(); i++) {
                double r = rdata.getValue(i);
                double g = gdata.getValue(i);
                double e = Math.exp(-sim.p2Truncated.u(r*r)/temperature);
                fw.write(String.format("%5.3f %22.15e %22.15e\n", r, g, e));
            }
            fw.close();
            
            if (computeP && computePMA) {
                fw = new FileWriter(functionsFile+"_mf.dat");
                rdata = meterF.getIndependentData(0);
                IData fdata = meterF.getData();
                Histogram f2hist = meterF.getHistogram2();
                double[] f2 = f2hist.getHistogram();
                for (int i=0; i<rdata.getLength(); i++) {
                    double r = rdata.getValue(i);
                    double mf = fdata.getValue(i);
                    if (Double.isNaN(mf)) continue;
                    double sdf = Math.sqrt(f2[i]-mf*mf);
                    double pf = -sim.p2Truncated.du(r*r)/r;
                    fw.write(String.format("%5.3f %22.15e %22.15e %22.15e\n", r, mf, sdf, pf));
                }
                fw.close();
            }
        }

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+(t2-t1)*0.001);
    }
    
    public static class LJMDParams extends ParameterBase {
        public int numAtoms = 100;
        public double temperature = 1.0;
        public double density = 0.01;
        public long numSteps = 1000000;
        public double rc = 4;
        public String functionsFile = null;
        public boolean computeP = true;
        public boolean computePMA = true;
        public boolean graphics = false;
        public int nBlocks = 100;
    }
}
