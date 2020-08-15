package etomica.osmoticvirial;

import etomica.action.BoxImposePbc;

import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPumpListener;
import etomica.data.histogram.HistogramSimple;
import etomica.data.meter.MeterNMolecules;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.DoubleRange;
import etomica.molecule.IMolecule;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.cell.PotentialMasterCellMixed;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 * Implements Ashton and Wilding method for calculation of osmotic virial coefficient (B2, B3)
 * for Hard-Sphere model as described in the paper DOI: 10.1063/1.4883718.
 */
public class AshtonWildingOsmoticVirial extends Simulation {

    protected Box box;
    protected Potential2 potential1, potential2, potential12;
    protected IntegratorMC integrator;
    protected MCMoveAtom mcMoveAtom;
    protected MCMoveInsertDelete mcMoveInsertDelete;
    protected MCMoveGeometricCluster mcMoveGeometricCluster;
    protected SpeciesSpheresMono species1, species2;
    

    /**
     * @param numAtoms no. of solute atoms in the box
     * @param vf reservoir volume fraction of solvent
     * @param q size of solvent divided by size of solute
     * @param computeIdeal whether to compute histograms for ideal gas
     */
    public AshtonWildingOsmoticVirial(int numAtoms, double vf, double q, boolean computeIdeal, double L, double GCfreq, boolean graphics){

        super(Space3D.getInstance());
//      setRandom(new RandomMersenneTwister(1));
//      PotentialMasterCell potentialMaster = new PotentialMasterCell(this, space);
        PotentialMaster potentialMaster;
        if(!computeIdeal) {
            potentialMaster = new PotentialMasterCellMixed(this, q);
        }
        else{
            potentialMaster = new PotentialMaster();
        }
        PotentialMasterCell pmc = potentialMaster instanceof PotentialMasterCell ? (PotentialMasterCell) potentialMaster : null;

        double sigma1 = 1.0; //solute
        double sigma2 = q * sigma1; //solvent
        double sigma12 = (sigma1+sigma2)/2;

        species1 = new SpeciesSpheresMono(this, space);
        species2 = new SpeciesSpheresMono(this, space);
        addSpecies(species1);
        addSpecies(species2);

        box = new Box(new BoundaryRectangularPeriodic(space, L * sigma1), space);
        addBox(box);
        box.setNMolecules(species1,numAtoms);

        integrator = new IntegratorMC(this, potentialMaster, box);
        this.getController2().addActivity(new ActivityIntegrate2(integrator));
        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        if (pmc != null) pmc.setCellRange(2);
        potentialMaster.setPotentialHard(true);

        System.out.println("vol: "+box.getBoundary().volume());

        if(computeIdeal) {
            potential1 = new P2Ideal(space);
            potential2 = new P2Ideal(space);
            potential12 = new P2Ideal(space);
            ((P2Ideal) potential1).setRange(1);
            System.out.println("P_ideal");
        }
        else{
            potential1 = new P2HardSphere(space, sigma1, false);
            potential2 = new P2HardSphere(space, sigma2, false);
//            potential2 = new P2Ideal(space);
            potential12 = new P2HardSphere(space, sigma12, false);
//            System.out.println("AO");
            if(vf != 0) {
                mcMoveInsertDelete = new MCMoveInsertDelete(potentialMaster, random, space);
                mcMoveInsertDelete.setSpecies(species2);
                double mu = (8 * vf - 9 * vf * vf + 3 * vf * vf * vf) / Math.pow((1 - vf), 3) + Math.log(6 * vf / (Math.PI * Math.pow(sigma2, 3))); //Configurational chemical potential from Carnahan–Starling equation of state
                System.out.println("mu " + mu + " muig " + Math.log(6 * vf / (Math.PI * Math.pow(sigma2, 3))));
//              mu = Double.MAX_VALUE;
                mcMoveInsertDelete.setMu(mu);
                integrator.getMoveManager().addMCMove(mcMoveInsertDelete);
                mcMoveGeometricCluster = new MCMoveGeometricCluster(pmc, space, random, integrator, species1);
                integrator.getMoveManager().addMCMove(mcMoveGeometricCluster);
            }
        }

//      potentialMaster.setCellRange(3);
//      potentialMaster.setRange(potential1.getRange());
//      potentialMaster.setPotentialHard(true);

        AtomType leafType1 = species1.getLeafType();
        AtomType leafType2 = species2.getLeafType();

        if (potentialMaster instanceof PotentialMasterCellMixed) {
            ((PotentialMasterCellMixed) potentialMaster).setHandledByCells(species2, true);
            ((PotentialMasterCellMixed) potentialMaster).addUnrangedPotential(leafType1, leafType1, potential1);
            ((PotentialMasterCellMixed) potentialMaster).addUnrangedPotential(leafType1, leafType2, potential12);
            potentialMaster.addPotential(potential2, new AtomType[]{leafType2, leafType2});
        } else {
            potentialMaster.addPotential(potential1, new AtomType[]{leafType1, leafType1});
            potentialMaster.addPotential(potential12, new AtomType[]{leafType1, leafType2});
            potentialMaster.addPotential(potential2, new AtomType[]{leafType2, leafType2});

            integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
            //To take care of positions updated by periodic boundary condition
        }

        if(!computeIdeal && vf != 0) {
            mcMoveGeometricCluster.setPotential(leafType1, leafType1, potential1);
            mcMoveGeometricCluster.setPotential(leafType1, leafType2, potential12);
            mcMoveGeometricCluster.setPotential(leafType2, leafType2, potential2);
            }


//      potentialMaster.addPotential(potential1, new AtomType[]{leafType1, leafType1});
//      potentialMaster.addPotential(potential12, new AtomType[]{leafType1, leafType2});
//      potentialMaster.addPotential(potential2, new AtomType[]{leafType2, leafType2});

        // add solvents to the boxes as though there are no solutes (and add solutes as though there
        // are no solvents).  then remove any solvents that overlap a solute
        int nSolvent = (int) (box.getBoundary().volume() * vf / (q * q * q * Math.PI / 6));
        box.setNMolecules(species2, nSolvent);

//      integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box, species1);
        configuration.initializeCoordinates(box, species2);

        if (!computeIdeal){
            final Set<IMolecule> overlaps = new HashSet<>();
            PotentialCalculation overlapCheck = new PotentialCalculation() {
                @Override
                public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
                    double u = potential.energy(atoms);
                    if (u < Double.POSITIVE_INFINITY) return;
                    IAtom a = atoms.get(0);
                    if (a.getType().getSpecies() == species1) a = atoms.get(1);
                    overlaps.add(a.getParentGroup());
                }
            };
            IteratorDirective id = new IteratorDirective();
            potentialMaster.calculate(box, id, overlapCheck);
            for (IMolecule m : overlaps) box.removeMolecule(m);
        }
        // GCfreq yields a simulation that spends roughly equal
        // CPU time doing GC moves and Displacement/Insert-Delete moves.  Set GCfreq to achieve desired result.

        if(!computeIdeal && vf != 0){
            integrator.getMoveManager().setFrequency(mcMoveGeometricCluster,2*GCfreq);
        }

        if (pmc != null){
            // we have cell listing.  initialize all that.
            integrator.getMoveEventManager().addListener(pmc.getNbrCellManager(box).makeMCMoveListener());
            pmc.getNbrCellManager(box).assignCellAll();
        }

//        integrator.setBox(box);
//        potentialMaster.getNbrCellManager(box).assignCellAll();
    }

    public static void main(String[] args) throws IOException {
        simParams params = new simParams();

        if(args.length > 0){
            ParseArgs.doParseArgs(params, args);
        }
        else{
            params.numAtoms = 2;
            params.numSteps = 1000000;
            params.nBlocks = 100;
            params.vf = 0.05;
            params.computeIdeal = false;
            params.sizeRatio = 0.33333333333333333333;
            params.L = 3;
            params.GCfreq = 3;
            params.graphics = false;
        }
        int numAtoms = params.numAtoms;
        long numSteps = params.numSteps;
        int nBlocks = params.nBlocks;
        double vf = params.vf;
        double q = params.sizeRatio;
        boolean computeIdeal = params.computeIdeal;
        double L = params.L;
        double GCfreq = params.GCfreq;
        boolean graphics = params.graphics;
        AccumulatorAverageFixed accNm = null;

        long numSamples = numSteps / numAtoms ;
        long samplesPerBlock = numSamples / nBlocks;
        if(samplesPerBlock == 0) samplesPerBlock = 1;

        System.out.println(numAtoms + " Atoms, "+ numSteps + " Steps" );
        System.out.println("Volume fraction: "+ vf);
        System.out.println("Size ratio: "+ q);
        System.out.println("nBlocks "+ nBlocks);
        System.out.println("system size: "+L+" sigma" );
        System.out.println("GC factor " + GCfreq);

        long t1 = System.currentTimeMillis();

        AshtonWildingOsmoticVirial sim = new AshtonWildingOsmoticVirial(numAtoms, vf, q, computeIdeal, L, GCfreq, graphics);

        if(graphics){
            final String appName = "Ashton-Wilding";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, appName, 3);

            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species1.getLeafType(), 1);
            ((DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species2.getLeafType(), q);
            simGraphic.makeAndDisplayFrame(appName);
            return;
        }
        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps/10);

sim.integrator.getMoveManager().setEquilibrating(false);

        MeterRminSpecies meterRmin = new MeterRminSpecies(sim.space, sim.box, sim.species1);
        AccumulatorHistogram accRmin = new AccumulatorHistogram(new HistogramSimple(new DoubleRange(0, L*sim.potential1.getRange())));
        DataPumpListener pumpRmin = new DataPumpListener(meterRmin,accRmin,numAtoms);
        sim.integrator.getEventManager().addListener(pumpRmin);

        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.species2);
        meterNMolecules.setBox(sim.box);
        accNm = new AccumulatorAverageFixed(samplesPerBlock);
        DataPumpListener pumpNm = new DataPumpListener(meterNMolecules, accNm);
        sim.integrator.getEventManager().addListener(pumpNm);
sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps);

        double[] histRmin = accRmin.getHistograms().getHistogram();
        double[] r = accRmin.getHistograms().xValues();
        System.out.println("\nWriting probabilities to test.txt file\n");
        FileWriter writer = new FileWriter("test.txt");
        for (int i = 0; i < histRmin.length; i++) {
            writer.write(r[i]+" "+histRmin[i]+"\n");
        }
        writer.close();

        double NmAvg = accNm.getData(AccumulatorAverage.AVERAGE).getValue(0);
        double NmErr = accNm.getData(AccumulatorAverage.ERROR).getValue(0);
        double NmCor = accNm.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);

        double VfBox = NmAvg*q*q*q*Math.PI/6/sim.box.getBoundary().volume();
        double VfErr = NmErr*q*q*q*Math.PI/6/sim.box.getBoundary().volume();

        System.out.print(String.format("No. of molecules avg: %13.6e  err: %11.4e   cor: % 4.2f\n", NmAvg, NmErr, NmCor));
        System.out.print(String.format("Volume frac avg: %13.6e  err: %11.4e\n", VfBox, VfErr));

        long t2 = System.currentTimeMillis();
        System.out.println("time: "+ (t2-t1)*0.001);

    }

    public static class simParams extends ParameterBase{
        public int numAtoms = 2;
        public long numSteps = 10000000;
        public int nBlocks = 100;
        public double vf = 0.0;
        public double sizeRatio = 0.33333333333333333333;
        public boolean computeIdeal = false;
        public double L = 3.0;
        public double GCfreq = 3;
        public boolean graphics = false;

    }
}
