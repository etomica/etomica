/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedvirial;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterRDFPC;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Histogram;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

public class MappedU extends Simulation {
    
    public SpeciesSpheresMono species;
    public IBox box;
    public IntegratorMC integrator;
    public MCMoveAtom move;
    public ActivityIntegrate activityIntegrate;
    public P2SoftSphericalTruncatedShifted p2TruncatedShifted;
    
    public MappedU(ISpace _space, int numAtoms, double temperature, double density, double rc) {
        super(_space);

	    //species and potentials
	    species = new SpeciesSpheresMono(this, space);
        addSpecies(species);
        
        //construct box
	    box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        PotentialMasterCell potentialMaster = new PotentialMasterCell(this, rc, space);
        potentialMaster.lrcMaster().setEnabled(false);
        
        //controller and integrator
	    integrator = new IntegratorMC(potentialMaster, random, temperature);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        move = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(move);

	    P2LennardJones potential = new P2LennardJones(space);
        p2TruncatedShifted = new P2SoftSphericalTruncatedShifted(space, potential, rc);
	    potentialMaster.addPotential(p2TruncatedShifted, new IAtomType[]{species.getLeafType(), species.getLeafType()});
	    
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        integrator.setBox(box);
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
            params.density = 1.0;
            params.numSteps = 10000000;
            params.rc = 3;
            params.numAtoms = 400;
            params.functionsFile = "0.10";
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

        System.out.println("Virial mapped average");
        System.out.println(numAtoms+" atoms, "+numSteps+" steps");
        System.out.println("density: "+density);
        System.out.println("temperature: "+temperature);
        System.out.println("cutoff: "+rc);
        System.out.println(nBlocks+" blocks");
        
        FileWriter fwr = new FileWriter("/usr/users/aksharag/workspace/Apps/Ushift/mapped_shifted.dat",true);

        ISpace space = Space.getInstance(3);

        MappedU sim = new MappedU(space, numAtoms, temperature, density, rc);
        
        if (graphics) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
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
           // rdfPlot.getPlot().setTitle("Radial Distribution Function");
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
 
        final MeterMappedU meterMappedU = new MeterMappedU(space, sim.integrator.getPotentialMaster(), sim.box, nBins);
        meterMappedU.getPotentialCalculation().setTemperature(sim.integrator.getTemperature(), sim.p2TruncatedShifted);
        final AccumulatorAverageFixed accMappedVirial = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpMappedU = new DataPumpListener(meterMappedU, accMappedVirial, numAtoms);
        if (computeUMA) sim.integrator.getEventManager().addListener(pumpMappedU);
        System.out.println("x0: "+meterMappedU.getPotentialCalculation().getX0());
        
        double qp_q = meterMappedU.getPotentialCalculation().getQP_Q();
        System.out.println("qp/q: "+qp_q);
        System.out.println("rho*db2/dbeta: "+qp_q*numAtoms/2);
        double qp = meterMappedU.getPotentialCalculation().getqp();
        System.out.println("qp: "+qp); 
        
        final MeterPotentialEnergy meterU = new MeterPotentialEnergy(sim.integrator.getPotentialMaster());
        meterU.setBox(sim.box);
        final AccumulatorAverageFixed accU = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpU = new DataPumpListener(meterU, accU, numAtoms);
        if (computeU) sim.integrator.getEventManager().addListener(pumpU);

        MeterMeanForce meterF = null;
        MeterRDFPC meterRDF = null;
        if (functionsFile != null) {
            int nbins = (int)Math.round(rc/0.01);
            meterF = new MeterMeanForce(space, sim.integrator.getPotentialMaster(), sim.p2TruncatedShifted, sim.box, nbins);
            if (computeU && computeUMA) sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterF, numAtoms));
    
            meterRDF = new MeterRDFPC(space, sim.integrator.getPotentialMaster(), sim.box);
            meterRDF.setBox(sim.box);
            meterRDF.getXDataSource().setNValues(nbins);
            meterRDF.getXDataSource().setXMax(rc);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, numAtoms));
        }

        sim.getController().actionPerformed();
        
        IData mappedAvg = accMappedVirial.getData(accMappedVirial.AVERAGE);
        IData mappedErr = accMappedVirial.getData(accMappedVirial.ERROR);
        IData mappedCor = accMappedVirial.getData(accMappedVirial.BLOCK_CORRELATION);
        double avg = mappedAvg.getValue(0);
        double err = mappedErr.getValue(0);
        double cor = mappedCor.getValue(0);
        
        double UavgMInt = (-1*qp_q*numAtoms/2)-(avg/numAtoms);
        
        if (computeUMA) {
        	System.out.print(String.format("avg: %13.6e   err: %11.4e   cor: % 4.2f\n", avg, err, cor));
        	System.out.println("Uavg extensive: "+ UavgMInt*numAtoms);
        	System.out.println("Uavg intensive: "+ UavgMInt);
        	System.out.println("error intensive: "+ err/numAtoms);
        }
             
        double UAvg = accU.getData(accU.AVERAGE).getValue(0);
        double UErr = accU.getData(accU.ERROR).getValue(0);
        double UCor = accU.getData(accU.BLOCK_CORRELATION).getValue(0);
        
        if (computeU){
        	System.out.print(String.format("Potential ext avg: %13.6e  err: %11.4e   cor: % 4.2f\n", UAvg, UErr, UCor));
        	System.out.print(String.format("Potential int avg:%13.6e  err int: %11.4e",UAvg/numAtoms,UErr/numAtoms));
        }
        
        fwr.append(numAtoms+" "+rc+" "+density+" "+temperature+" "+UavgMInt+" "+UAvg/numAtoms+" "+err/numAtoms+" "+UErr/numAtoms+"\n");
        
        if (functionsFile != null) {
            FileWriter fw = new FileWriter(functionsFile+"_gr.dat");
            IData rdata = meterRDF.getIndependentData(0);
            IData gdata = meterRDF.getData();
            for (int i=0; i<rdata.getLength(); i++) {
                double r = rdata.getValue(i);
                double g = gdata.getValue(i);
                double e = Math.exp(-sim.p2TruncatedShifted.u(r*r)/temperature);
                fw.write(String.format("%5.3f %22.15e %22.15e\n", r, g, e));
            }
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
                    double pf = -sim.p2TruncatedShifted.du(r*r)/r;
                    fw.write(String.format("%5.3f %22.15e %22.15e %22.15e\n", r, mf, sdf, pf));
                }
                fw.close();
            }
        }

        fwr.close();
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
        public int nBlocks = 100;
    }
}
