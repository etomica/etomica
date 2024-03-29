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
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterRDF;
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

public class MappedU extends Simulation {

    public SpeciesGeneral species;
    public Box box;
    public IntegratorMC integrator;
    public MCMoveAtom move;

    public P2SoftSphericalTruncated p2Truncated;

    public MappedU(Space _space, int numAtoms, double temperature, double density, double rc) {
        super(_space);

        //species and potentials
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        //construct box
        box = this.makeBox();
        box.setNMolecules(species, numAtoms);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        PotentialMasterCell potentialMaster = new PotentialMasterCell(getSpeciesManager(), box, 2, BondingInfo.noBonding());
        potentialMaster.doAllTruncationCorrection = false;

        //controller and integrator
        integrator = new IntegratorMC(potentialMaster, random, temperature, box);
        this.getController().addActivity(new ActivityIntegrate(integrator));
        move = new MCMoveAtom(random, potentialMaster, box);
        integrator.getMoveManager().addMCMove(move);

        P2LennardJones potential = new P2LennardJones();
        p2Truncated = new P2SoftSphericalTruncated(potential, rc);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), p2Truncated);

        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
    }

    public static void main(String[] args) throws IOException {

        LJMDParams params = new LJMDParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.temperature = 2.0;
            params.density = 0.01;
            params.numSteps = 10000000;
            params.rc = 2.5;
            params.numAtoms = 1000;
            params.functionsFile = null;
        }

        int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double density = params.density;
        long numSteps = params.numSteps;
        double rc = params.rc;
        String functionsFile = params.functionsFile;
        boolean computeU = params.computeU;
        boolean computeUMA = params.computeUMA;
        boolean graphics = params.graphics;
        int nBlocks = params.nBlocks;
        double halfBoxlength = 0;
        double qp_q = 0;
        AccumulatorAverageFixed accMappedVirial = null;
        AccumulatorAverageFixed accU = null;

        System.out.println("Virial mapped average Potential");
        System.out.println(numAtoms+" atoms, "+numSteps+" steps");
        System.out.println("density: "+density);
        System.out.println("temperature: "+temperature);
        System.out.println("cutoff: "+rc);
        System.out.println(nBlocks+" blocks");
         System.out.println("vshift: -p2.u(vCut*vCut)+0.0494908");

        // FileWriter fwr = new FileWriter("/usr/users/aksharag/workspace/Apps/Ushift/mapped_shifted.dat",true);

        Space space = Space.getInstance(3);

        MappedU sim = new MappedU(space, numAtoms, temperature, density, rc);
        halfBoxlength = Math.cbrt(numAtoms/density) / 2;

        System.out.println("half box length "+halfBoxlength);

        if (graphics) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.makeAndDisplayFrame();

            ArrayList<DataPump> dataStreamPumps = simGraphic.getController().getDataStreamPumps();

            MeterRDF meterRDF = new MeterRDF(space);
            meterRDF.setBox(sim.box);
            int nbins = (int)Math.round(rc/0.01);
            meterRDF.getXDataSource().setNValues(nbins);
            meterRDF.getXDataSource().setXMax(halfBoxlength);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, numAtoms));

            DisplayPlot rdfPlot = new DisplayPlot();
            DataPump rdfPump = new DataPump(meterRDF,rdfPlot.getDataSet().makeDataSink());
            IntegratorListenerAction rdfPumpListener = new IntegratorListenerAction(rdfPump);
            sim.integrator.getEventManager().addListener(rdfPumpListener);
            rdfPumpListener.setInterval(10);
            dataStreamPumps.add(rdfPump);

            rdfPlot.setDoLegend(false);
            // rdfPlot.getPlot().setTitle("Radial Distribution Function");
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

        if (computeUMA){
            final MeterMappedU meterMappedU = new MeterMappedU(sim.integrator.getPotentialCompute(), sim.box, nBins);
            meterMappedU.getPotentialCallback().setTemperature(sim.integrator.getTemperature(), sim.p2Truncated);
            accMappedVirial = new AccumulatorAverageFixed(samplesPerBlock);
            DataPumpListener pumpMappedU = new DataPumpListener(meterMappedU, accMappedVirial, numAtoms);
            if (computeUMA) sim.integrator.getEventManager().addListener(pumpMappedU);
            System.out.println("x0: "+meterMappedU.getPotentialCallback().getX0());

            qp_q = meterMappedU.getPotentialCallback().getQP_Q();
            System.out.println("qp/q: "+qp_q);
            System.out.println("rho*db2/dbeta: "+qp_q*numAtoms/2);
            double qp = meterMappedU.getPotentialCallback().getqp();
            System.out.println("qp: "+qp);
        }


        if (computeU){
            final MeterPotentialEnergyFromIntegrator meterU = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            accU = new AccumulatorAverageFixed(samplesPerBlock);
            DataPumpListener pumpU = new DataPumpListener(meterU, accU, 20);
            sim.integrator.getEventManager().addListener(pumpU);
        }

        MeterMeanForce meterF = null;
        MeterRDF meterRDF = null;
        if (functionsFile != null) {
            int nbins = (int)Math.floor(halfBoxlength/0.01);
            double eqncutoff = nbins * 0.01;
            meterF = new MeterMeanForce(sim.integrator.getPotentialCompute(), sim.p2Truncated, sim.box, nbins);
            if (computeU && computeUMA) sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterF, numAtoms));

            meterRDF = new MeterRDF(space);
            meterRDF.setBox(sim.box);
            meterRDF.getXDataSource().setNValues(nbins);
            meterRDF.getXDataSource().setXMax(eqncutoff);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, numAtoms));
        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));



        if (computeUMA) {

            IData mappedAvg = accMappedVirial.getData(accMappedVirial.AVERAGE);
            IData mappedErr = accMappedVirial.getData(accMappedVirial.ERROR);
            IData mappedCor = accMappedVirial.getData(accMappedVirial.BLOCK_CORRELATION);

            double avg = mappedAvg.getValue(0);
            double err = mappedErr.getValue(0);
            double cor = mappedCor.getValue(0);

            double UavgMInt = (-1*qp_q*numAtoms/2)-(avg/numAtoms);

            System.out.printf("avg: %13.6e   err: %11.4e   cor: % 4.2f\n", avg, err, cor);
            System.out.println("UavgM extensive: "+ UavgMInt*numAtoms);
            System.out.println("UavgM intensive: "+ UavgMInt);
            System.out.println("errorM intensive: "+ err/numAtoms);
        }



        if (computeU){

            double UAvg = accU.getData(accU.AVERAGE).getValue(0);
            double UErr = accU.getData(accU.ERROR).getValue(0);
            double UCor = accU.getData(accU.BLOCK_CORRELATION).getValue(0);

            System.out.println("Uavg intensive "+UAvg/numAtoms);
            System.out.println("err intensive "+UErr/numAtoms);
            System.out.printf("Potential ext avg: %13.6e  err: %11.4e   cor: % 4.2f\n", UAvg, UErr, UCor);
        }

        //  fwr.append(numAtoms+" "+rc+" "+density+" "+temperature+" "+UavgMInt+" "+UAvg/numAtoms+" "+err/numAtoms+" "+UErr/numAtoms+"\n");

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
            //System.out.println(rdata.getLength());
            fw.close();

            if (computeU && computeUMA) {
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

        //   fwr.close();
        long t2 = System.currentTimeMillis();
        System.out.println("\n time: "+(t2-t1)*0.001);

    }

    public static class LJMDParams extends ParameterBase {
        public int numAtoms = 100;
        public double temperature = 1.0;
        public double density = 0.01;
        public long numSteps = 1000000;
        public double rc = 4;
        public String functionsFile = null;
        public boolean computeU = true;
        public boolean computeUMA = true;
        public boolean graphics = false;
        public int nBlocks = 1000;
    }
}
