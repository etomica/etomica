
package etomica.mappedRdf2;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterRDF;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.io.IOException;

public class MappedRdf extends Simulation {

    public SpeciesSpheresMono species;
    public Box box;
    public IntegratorMC integrator;
    public MCMoveAtom move;
    public ActivityIntegrate activityIntegrate;
    public P2SoftSphericalTruncated p2Truncated;


    public MappedRdf(Space _space, int numAtoms, double temperature, double density, double rc) {
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
        p2Truncated = new P2SoftSphericalTruncated(space, potential, rc);
        potentialMaster.addPotential(p2Truncated, new IAtomType[]{species.getLeafType(), species.getLeafType()});

        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        integrator.setBox(box);
        potentialMaster.setCellRange(2);

        potentialMaster.getNbrCellManager(box).assignCellAll();

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }

    public static void main(String[] args) throws IOException {

        MappedRdf.LJMDParams params = new MappedRdf.LJMDParams();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.temperature = 5000.0;
            params.density = 0.01;
            params.numSteps = 300000;
            params.rc = 4;
            params.numAtoms = 500;

        }

        int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double density = params.density;
        long numSteps = params.numSteps;
        double rc = params.rc;
        boolean graphics = false;
        boolean computeR = params.computeR;
        boolean computeRMA = params.computeRMA;

        int nBlocks = params.nBlocks;


        Space space = Space.getInstance(3);

        MappedRdf sim = new MappedRdf(space, numAtoms, temperature, density, rc);

        MeterRDF meterRDF = null;
        MeterMappedRdf meterMappedRdf = null;

        double halfBoxlength = sim.box.getBoundary().getBoxSize().getX(0) / 2;
        int nbins = (int) Math.floor(rc / 0.01);
        double eqncutoff = nbins * 0.01;


        if (graphics) {
            meterRDF = new MeterRDF(space);
            meterRDF.setBox(sim.box);
            meterRDF.getXDataSource().setNValues(nbins);
            meterRDF.getXDataSource().setXMax(eqncutoff);

            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, numAtoms));

            DisplayPlot rdfPlot = new DisplayPlot();
            DataPump rdfPump = new DataPump(meterRDF, rdfPlot.getDataSet().makeDataSink());
            IntegratorListenerAction rdfPumpListener = new IntegratorListenerAction(rdfPump);
            sim.integrator.getEventManager().addListener(rdfPumpListener);
            rdfPumpListener.setInterval(10 * numAtoms);

            SimulationGraphic gsim = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());

            rdfPlot.setDoLegend(false);
            rdfPlot.getPlot().setTitle("Radial Distribution Function");
            rdfPlot.setLabel("RDF");

            gsim.add(rdfPlot);
            gsim.makeAndDisplayFrame();

            return;

        }

        sim.activityIntegrate.setMaxSteps(numSteps / 10);
        sim.activityIntegrate.actionPerformed();
        sim.activityIntegrate.setMaxSteps(numSteps);

        sim.integrator.getMoveManager().setEquilibrating(false);

        // if(computeR){

        meterRDF = new MeterRDF(space,true);
        meterRDF.setBox(sim.box);
        meterRDF.getXDataSource().setNValues(nbins);
        meterRDF.getXDataSource().setXMax(eqncutoff);

        meterMappedRdf = new MeterMappedRdf(space, sim.integrator.getPotentialMaster(), sim.box, nbins);
        meterMappedRdf.setBox(sim.box);
        meterMappedRdf.getXDataSource().setNValues(nbins);
        meterMappedRdf.getXDataSource().setXMax(eqncutoff);

        meterMappedRdf.getPotentialCalculation().setPotential(sim.p2Truncated);
        meterMappedRdf.getPotentialCalculation().setTemperature(temperature);

        AccumulatorAverageFixed accmap = new AccumulatorAverageFixed(numSteps/(nBlocks*numAtoms));
        DataPumpListener map = new DataPumpListener(meterMappedRdf,accmap,numAtoms);
        sim.integrator.getEventManager().addListener(map);

        AccumulatorAverageFixed acccon = new AccumulatorAverageFixed(numSteps/(nBlocks*numAtoms));
        DataPumpListener con = new DataPumpListener(meterRDF,acccon,numAtoms);
        sim.integrator.getEventManager().addListener(con);

        sim.activityIntegrate.actionPerformed();

        IData rdata = meterRDF.getIndependentData(0);
        IData gdata = acccon.getData(acccon.AVERAGE);
        IData gdataerr = acccon.getData(acccon.ERROR);
        IData rmdata = meterMappedRdf.getIndependentData(0);
        IData gmdata = accmap.getData(accmap.AVERAGE);
        IData gmdataerr = accmap.getData(accmap.ERROR);

        for (int i = 0; i < rdata.getLength(); i++)

        {
            double r = rdata.getValue(i);
            double g = gdata.getValue(i);
            double gerr = gdataerr.getValue(i);
            double gm = gmdata.getValue(i);
            double gmerr = gmdataerr.getValue(i);

            // double e = Math.exp(-sim.p2Truncated.u(r * r) / temperature);
            System.out.println(r + " " + g + " " + gerr + " "  + gm+  " " + gmerr);
        }

        // }

        //   if(computeRMA){


        // fw.close();

//        FileWriter fw = new FileWriter("rdfDiff.dat");
//        for (int i = 0; i < rmdata.getLength(); i++) {
//            double rm = rmdata.getValue(i);
//            double gdiff = gdata.getValue(i) - gmdata.getValue(i);
//            System.out.println("rm " + rm + " gDiff " + gdiff);
//            fw.write(rm + " " + gdiff + "\n");
//        }
//        fw.close();

    }

    // }


    public static class LJMDParams extends ParameterBase {
        public int numAtoms = 100;
        public double temperature = 1.0;
        public double density = 0.01;
        public long numSteps = 1000000;
        public double rc = 4;
        public int nBlocks = 100;
        public boolean computeR = false;
        public boolean computeRMA = true;
    }
}
