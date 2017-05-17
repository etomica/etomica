
package etomica.mappedRdf;

import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.IAtomType;
import etomica.data.IData;
import etomica.data.meter.MeterRDF;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.listener.IntegratorListenerAction;

import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.config.ConfigurationLattice;

import etomica.box.Box;
import etomica.action.BoxInflate;
import etomica.lattice.LatticeCubicFcc;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Created by aksharag on 5/15/17.
 */
public class MappedRdf extends Simulation {

    public SpeciesSpheresMono species;
    public IBox box;
    public IntegratorMC integrator;
    public MCMoveAtom move;
    public ActivityIntegrate activityIntegrate;
    public P2SoftSphericalTruncated p2Truncated;


    public MappedRdf(ISpace _space, int numAtoms, double temperature, double density, double rc) {
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
        }

        else {
            params.temperature = 2.0;
            params.density = 0.01;
            params.numSteps = 10000000;
            params.rc = 2.5;
            params.numAtoms = 1000;

        }

        int numAtoms = params.numAtoms;
        double temperature = params.temperature;
        double density = params.density;
        long numSteps = params.numSteps;
        double rc = params.rc;

        int nBlocks = params.nBlocks;
        double halfBoxlength = 0;

        ISpace space = Space.getInstance(3);

        MappedRdf sim = new MappedRdf(space, numAtoms, temperature, density, rc);
        MeterRDF meterRDF = null;

        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.integrator.getMoveManager().setEquilibrating(false);

        halfBoxlength = sim.box.getBoundary().getBoxSize().getX(0) /2;

        int nbins = (int) Math.floor(halfBoxlength / 0.01);
        double eqncutoff = nbins * 0.01;

        meterRDF = new MeterRDF(space);
        meterRDF.setBox(sim.box);
        meterRDF.getXDataSource().setNValues(nbins);
        meterRDF.getXDataSource().setXMax(eqncutoff);


        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterRDF, numAtoms));

       // sim.getController().actionPerformed();



        IData rdata = meterRDF.getIndependentData(0);
        IData gdata = meterRDF.getData();

        for (int i = 0; i < rdata.getLength(); i++)

        {
            double r = rdata.getValue(i);
            double g = gdata.getValue(i);
            double e = Math.exp(-sim.p2Truncated.u(r * r) / temperature);
            System.out.println("r "+r+" g "+g);
        }

    }
    public static class LJMDParams extends ParameterBase {
        public int numAtoms = 100;
        public double temperature = 1.0;
        public double density = 0.01;
        public long numSteps = 1000000;
        public double rc = 4;
        public int nBlocks = 1000;
    }
}

