/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.histogram.Histogram;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterRDFPC;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.integrator.IntegratorListenerAction;
import etomica.mappedvirial.PotentialCalculationMappedVirialV.VFunc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class MappedVirialLJVGr extends Simulation {
    
    public SpeciesSpheresMono species;
    public Box box;
    public IntegratorMC integrator;
    public MCMoveAtom move;
    public ActivityIntegrate activityIntegrate;
    public P2SoftSphericalTruncated p2Truncated;
    
    public MappedVirialLJVGr(Space _space, int numAtoms, double temperature, double density, double rc) {
        super(_space);
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, rc, space);
        potentialMaster.lrcMaster().setEnabled(false);

        //controller and integrator
        box = this.makeBox();
        integrator = new IntegratorMC(potentialMaster, random, temperature, box);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        move = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(move);

        //species and potentials
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        P2LennardJones potential = new P2LennardJones(space);
        p2Truncated = new P2SoftSphericalTruncated(space, potential, rc);
        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});

        //construct box
        box.setNMolecules(species, numAtoms);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        potentialMaster.setCellRange(2);

        potentialMaster.getNbrCellManager(box).assignCellAll();

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }
    
    public static void main(String[] args) throws IOException {
        
        LJMDParams params = new LJMDParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.temperature = 2.0;
            params.density = 0.30;
            params.numSteps = 1000000;
            params.rc = 4;
            params.numAtoms = 1000;
            params.functionsFile = "0.30";
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
        VSource vSource = params.vSource;

        System.out.println("Virial mapped average");
        System.out.println(numAtoms+" atoms, "+numSteps+" steps");
        System.out.println("density: "+density);
        System.out.println("temperature: "+temperature);
        System.out.println("cutoff: "+rc);
        System.out.println(nBlocks+" blocks");

        Space space = Space.getInstance(3);

        MappedVirialLJVGr sim = new MappedVirialLJVGr(space, numAtoms, temperature, density, rc);
        
        if (graphics) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.makeAndDisplayFrame();

            ArrayList<DataPump> dataStreamPumps = simGraphic.getController().getDataStreamPumps();
            
            MeterRDFPC meterRDF = new MeterRDFPC(space, sim.integrator.getPotentialMaster(), sim.box);
            meterRDF.setBox(sim.box);
            int nbins = (int)Math.round(rc/0.01);
            meterRDF.getXDataSource().setNValues(nbins);
            meterRDF.getXDataSource().setXMax(rc);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, numAtoms));

            DisplayPlot rdfPlot = new DisplayPlot();
            DataPump rdfPump = new DataPump(meterRDF,rdfPlot.getDataSet().makeDataSink());
            IntegratorListenerAction rdfPumpListener = new IntegratorListenerAction(rdfPump);
            sim.integrator.getEventManager().addListener(rdfPumpListener);
            rdfPumpListener.setInterval(10);
            dataStreamPumps.add(rdfPump);
            
            rdfPlot.setDoLegend(false);
            rdfPlot.getPlot().setTitle("Radial Distribution Function");
            rdfPlot.setLabel("RDF");
            simGraphic.add(rdfPlot);

            return;
        }
        
        long t1 = System.currentTimeMillis();
        
        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.integrator.getMoveManager().setEquilibrating(false);
        
        int nBins = 1000000;
        long numSamples = numSteps/numAtoms;
        long samplesPerBlock = numSamples/nBlocks;
        if (samplesPerBlock == 0) samplesPerBlock = 1;

        final MeterMappedVirialV meterMappedVirial = new MeterMappedVirialV(space, sim.integrator.getPotentialMaster(), sim.box, nBins);
        meterMappedVirial.getPotentialCalculation().setTemperature(sim.integrator.getTemperature(), sim.p2Truncated);
        if (vSource == VSource.SINE) {
            VFunc v = new PotentialCalculationMappedVirialV.USine(sim.p2Truncated, new double[][]{{0.87624,2.17775},{0.684801,0.647015},{-0.748087,1.1296},{0.692743,-1.07915}}, temperature);
            meterMappedVirial.getPotentialCalculation().setVFunc(v);
        }
        final AccumulatorAverageFixed accMappedVirial = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpMappedVirial = new DataPumpListener(meterMappedVirial, accMappedVirial, numAtoms);
        if (computePMA) sim.integrator.getEventManager().addListener(pumpMappedVirial);
        System.out.println("x0: "+meterMappedVirial.getPotentialCalculation().getX0());
        double qu = meterMappedVirial.getPotentialCalculation().getQU();
        System.out.println("qu: "+qu);
        double q = meterMappedVirial.getPotentialCalculation().getQ();
        System.out.println("q: "+q);

        final MeterPressure meterP = new MeterPressure(space);
        meterP.setIntegrator(sim.integrator);
        final AccumulatorAverageFixed accP = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpP = new DataPumpListener(meterP, accP, numAtoms);
        if (computeP) sim.integrator.getEventManager().addListener(pumpP);

        MeterMeanForce meterF = null;
        MeterRDFPC meterRDF = null;
        if (functionsFile != null) {
            int nbins = (int)Math.round(rc/0.01);
            meterF = new MeterMeanForce(space, sim.integrator.getPotentialMaster(), sim.p2Truncated, sim.box, nbins);
            if (computeP && computePMA) sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterF, numAtoms));
    
            meterRDF = new MeterRDFPC(space, sim.integrator.getPotentialMaster(), sim.box);
            meterRDF.setBox(sim.box);
            meterRDF.getXDataSource().setNValues(nbins);
            meterRDF.getXDataSource().setXMax(rc);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, numAtoms));
        }

        sim.getController().actionPerformed();

        IData mappedAvg = accMappedVirial.getData(AccumulatorAverage.AVERAGE);
        IData mappedErr = accMappedVirial.getData(AccumulatorAverage.ERROR);
        IData mappedCor = accMappedVirial.getData(AccumulatorAverage.BLOCK_CORRELATION);
        double avg = mappedAvg.getValue(0);
        double err = mappedErr.getValue(0);
        double cor = mappedCor.getValue(0);
        if (computePMA) System.out.print(String.format("avg: %13.6e   err: %11.4e   cor: % 4.2f\n", avg, err, cor));

        double pAvg = accP.getData(AccumulatorAverage.AVERAGE).getValue(0);
        double pErr = accP.getData(AccumulatorAverage.ERROR).getValue(0);
        double pCor = accP.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);
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

    enum VSource {U_SHIFT, TANH, SINE}

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
        public VSource vSource = VSource.U_SHIFT; 
    }
}
