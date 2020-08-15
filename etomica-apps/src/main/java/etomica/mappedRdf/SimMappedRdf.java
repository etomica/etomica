package etomica.mappedRdf;

import etomica.action.BoxInflate;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterRDF;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Calculates pair distribution using histograms and mapped averaging
 */
public class SimMappedRdf extends Simulation {

    public SpeciesSpheresMono species;
    public Box box;
    public IntegratorMC integrator;
    public MCMoveAtom move;

    //   public P2SoftSphericalTruncatedForceShifted p2Truncated;
    public P2SoftSphericalTruncatedShifted p2Truncated;

    public SimMappedRdf(Space _space, int numAtoms, double temperature, double density, double rc) {
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
        integrator = new IntegratorMC(potentialMaster, random, temperature, box);
        this.getController().addActivity(new ActivityIntegrate(integrator));
        move = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(move);

        P2LennardJones potential = new P2LennardJones(space);
        //   p2Truncated = new P2SoftSphericalTruncatedForceShifted(space, potential, rc);
        p2Truncated = new P2SoftSphericalTruncatedShifted(space, potential, rc);

        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});

        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        potentialMaster.setCellRange(2);

        potentialMaster.getNbrCellManager(box).assignCellAll();

        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }

    public static void main(String[] args) {

        LJMDParams params = new LJMDParams();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.numSteps = 10000;
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

        SimMappedRdf sim = new SimMappedRdf(space, numAtoms, temperature, density, rc);

        MeterRDF  meterRDF = null;
        MeterMappedRdf meterMappedRdf = null;

        double halfBoxlength = sim.box.getBoundary().getBoxSize().getX(0) / 2;
        int nbins = 100;
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

            SimulationGraphic gsim = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

            rdfPlot.setDoLegend(false);
            rdfPlot.getPlot().setTitle("Radial Distribution Function");
            rdfPlot.setLabel("RDF");

            gsim.add(rdfPlot);
            gsim.makeAndDisplayFrame();

            return;

        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), numSteps / 10);

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

        AccumulatorAverageFixed accmap = new AccumulatorAverageFixed(numSteps/(nBlocks*numAtoms));
        DataPumpListener map = new DataPumpListener(meterMappedRdf,accmap,numAtoms);
        sim.integrator.getEventManager().addListener(map);

        AccumulatorAverageFixed acccon = new AccumulatorAverageFixed(numSteps/(nBlocks*numAtoms));
        DataPumpListener con = new DataPumpListener(meterRDF,acccon,numAtoms);
        sim.integrator.getEventManager().addListener(con);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), numSteps);

        IData rdata = meterRDF.getIndependentData(0);
        IData gdata = acccon.getData(acccon.AVERAGE);
        IData gdataerr = acccon.getData(acccon.ERROR);

        //  IData rmdata = meterMappedRdf.getIndependentData(0);
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
            System.out.println(r + " " + g*params.numAtoms*(params.numAtoms-1)/((params.numAtoms/params.density)*(params.numAtoms/params.density)) + " " + gerr*params.numAtoms*(params.numAtoms-1)/((params.numAtoms/params.density)*(params.numAtoms/params.density)) + " "  + gm+  " " + gmerr);
        }

    }

    public static class LJMDParams extends ParameterBase {
        public int numAtoms = 100;
        public double temperature = 5.0;
        public double density = 0.125;
        public long numSteps =  100;
        public double rc = 3;
        public int nBlocks = 100;
        public boolean computeR = false;
        public boolean computeRMA = true;
        public double rcforHandfinmap = 3;

    }
}
