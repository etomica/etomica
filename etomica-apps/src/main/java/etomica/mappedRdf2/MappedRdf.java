
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
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.P2SoftSphericalTruncatedShifted;
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
 //   public P2SoftSphericalTruncatedForceShifted p2Truncated;
 public P2SoftSphericalTruncatedShifted p2Truncated;

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
     //   p2Truncated = new P2SoftSphericalTruncatedForceShifted(space, potential, rc);
        p2Truncated = new P2SoftSphericalTruncatedShifted(space, potential, rc);

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
            params.numSteps = 100000000;
        }

        int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double density = params.density;
        long numSteps = params.numSteps;
        double rc = params.rc;
        double rcforHandfinmap = params.rcforHandfinmap;
        boolean graphics = false;
        boolean computeR = params.computeR;
        boolean computeRMA = params.computeRMA;

        int nBlocks = params.nBlocks;


        Space space = Space.getInstance(3);

        MappedRdf sim = new MappedRdf(space, numAtoms, temperature, density, rc);

        MeterRDF  meterRDF = null;
        MeterMappedRdf meterMappedRdf = null;
        MeterMappedRdf2 meterMappedRdf2 = null;
        MeterMappedRdf3 meterMappedRdf3 = null;

        double halfBoxlength = sim.box.getBoundary().getBoxSize().getX(0) / 2;
        int nbins = 1000;
        double eqncutoff = halfBoxlength;


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

        meterMappedRdf = new MeterMappedRdf(params.rcforHandfinmap,space, sim.integrator.getPotentialMaster(), sim.box, nbins,params.density);
        meterMappedRdf.setBox(sim.box);
        meterMappedRdf.getXDataSource().setNValues(nbins);
        meterMappedRdf.getXDataSource().setXMax(eqncutoff);
        meterMappedRdf.getPotentialCalculation().setPotential(sim.p2Truncated);
        meterMappedRdf.getPotentialCalculation().setTemperature(temperature);

        meterMappedRdf2 = new MeterMappedRdf2(params.rcforHandfinmap,space, sim.integrator.getPotentialMaster(), sim.box, nbins,params.density);
        meterMappedRdf2.setBox(sim.box);
        meterMappedRdf2.getXDataSource().setNValues(nbins);
        meterMappedRdf2.getXDataSource().setXMax(eqncutoff);
        meterMappedRdf2.getPotentialCalculation().setPotential(sim.p2Truncated);
        meterMappedRdf2.getPotentialCalculation().setTemperature(temperature);

        meterMappedRdf3 = new MeterMappedRdf3(params.rcforHandfinmap,space, sim.integrator.getPotentialMaster(), sim.box, nbins,params.density);
        meterMappedRdf3.setBox(sim.box);
        meterMappedRdf3.getXDataSource().setNValues(nbins);
        meterMappedRdf3.getXDataSource().setXMax(eqncutoff);
        meterMappedRdf3.getPotentialCalculation().setPotential(sim.p2Truncated);
        meterMappedRdf3.getPotentialCalculation().setTemperature(temperature);


        AccumulatorAverageFixed accmap = new AccumulatorAverageFixed(numSteps/(nBlocks*numAtoms));
        DataPumpListener map = new DataPumpListener(meterMappedRdf,accmap,numAtoms);
        sim.integrator.getEventManager().addListener(map);

        AccumulatorAverageFixed acccon = new AccumulatorAverageFixed(numSteps/(nBlocks*numAtoms));
        DataPumpListener con = new DataPumpListener(meterRDF,acccon,numAtoms);
        sim.integrator.getEventManager().addListener(con);

        AccumulatorAverageFixed accwithout3rdterm = new AccumulatorAverageFixed(numSteps/(nBlocks*numAtoms));
        DataPumpListener without3rdterm = new DataPumpListener(meterMappedRdf2,accwithout3rdterm,numAtoms);
        sim.integrator.getEventManager().addListener(without3rdterm);


        AccumulatorAverageFixed acc3rdterm = new AccumulatorAverageFixed(numSteps/(nBlocks*numAtoms));
        DataPumpListener term3rd = new DataPumpListener(meterMappedRdf3,acc3rdterm,numAtoms);
        sim.integrator.getEventManager().addListener(term3rd);
        sim.activityIntegrate.actionPerformed();

        IData rdata = meterRDF.getIndependentData(0);
        IData gdata = acccon.getData(acccon.AVERAGE);
        IData gdataerr = acccon.getData(acccon.ERROR);
        IData gwithout3rdtermdata = accwithout3rdterm.getData(accwithout3rdterm.AVERAGE);
        IData gwithout3rdtermdataerr = accwithout3rdterm.getData(accwithout3rdterm.ERROR);
      //  IData rmdata = meterMappedRdf.getIndependentData(0);
        IData gmdata = accmap.getData(accmap.AVERAGE);
        IData gmdataerr = accmap.getData(accmap.ERROR);
        IData g3rdtermdata = acc3rdterm.getData(acc3rdterm.AVERAGE);
        IData g3rdtermdataerr = acc3rdterm.getData(acc3rdterm.ERROR);

        for (int i = 0; i < rdata.getLength(); i++)

        {
            double r = rdata.getValue(i);
            double g = gdata.getValue(i);
            double gerr = gdataerr.getValue(i);
            double gwithout3rdterm = gwithout3rdtermdata.getValue(i);
            double gwithout3rdtermerr = gwithout3rdtermdataerr.getValue(i);
            double gm = gmdata.getValue(i);
            double gmerr = gmdataerr.getValue(i);
            double g3rdterm = g3rdtermdata.getValue(i);
            double g3rdtermerr = g3rdtermdataerr.getValue(i);

            // double e = Math.exp(-sim.p2Truncated.u(r * r) / temperature);
            System.out.println(r + " " + g*params.numAtoms*(params.numAtoms-1)/((params.numAtoms/params.density)*(params.numAtoms/params.density)) + " " + gerr*params.numAtoms*(params.numAtoms-1)/((params.numAtoms/params.density)*(params.numAtoms/params.density)) + " "  + gm+  " " + gmerr+ " "+ gwithout3rdterm +  " " + gwithout3rdtermerr+ " "+ g3rdterm +  " " + g3rdtermerr);
        }

    }

    public static class LJMDParams extends ParameterBase {
        public int numAtoms = 100;
        public double temperature = 5.0;
        public double density = 1;
        public long numSteps =  10000;
        public double rc = 2;
        public int nBlocks = 100;
        public boolean computeR = false;
        public boolean computeRMA = true;
        public double rcforHandfinmap = 2;

    }
}
