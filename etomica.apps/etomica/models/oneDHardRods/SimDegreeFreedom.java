package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPump;
import etomica.data.DataSplitter;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.MCMoveAtomCoupled;
import etomica.normalmode.NormalModes;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.P2XOrder;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.DoubleRange;

/**
 * MD simulation of hard spheres in 1D or 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 */
public class SimDegreeFreedom extends Simulation {

    public SimDegreeFreedom(Space _space, int numAtoms, double density) {
        super(_space, true);
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(this, potentialMaster);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        Potential potential = new P2HardSphere(space, 1.0, false);
        potential = new P2XOrder(space, (Potential2HardSpherical)potential);
        AtomTypeSphere sphereType = (AtomTypeSphere)species.getLeafType();
        potentialMaster.addPotential(potential, new IAtomType[] { sphereType,
                sphereType });

        mcmove = new MCMoveAtomCoupled(potentialMaster, random, space);
        mcmove.setPotential((Potential2)potential);
        mcmove.setBox(box);
        integrator.getMoveManager().addMCMove(mcmove);
        
        int nCells;
        Basis basis;
        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            nCells = numAtoms;
            bdry = new BoundaryRectangularPeriodic(space, numAtoms/density);
//            ((IntegratorHard) integrator).setNullPotential(new P1HardPeriodic(space), sphereType);
            basis = new BasisMonatomic(space);
        } else {
            primitive = new PrimitiveCubic(space, 1);
            double v = primitive.unitCell().getVolume();
            primitive.scaleSize(Math.pow(v*density/4,-1.0/3.0));
            nCells = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            bdry = new BoundaryDeformableLattice(primitive, new int[]{nCells,nCells,nCells});
            basis = new BasisCubicFcc();
        }
        box.setBoundary(bdry);
        integrator.setBox(box);
        
        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{nCells});
        int coordinateDim = coordinateDefinition.getCoordinateDim();
        
        if (space.D() == 1) {
            nm = new NormalModes1DHR(space.D());
        } else if (space.D() == 2) {
            waveVectorFactory = null;
        } else {
            waveVectorFactory = new WaveVectorFactorySimple(primitive, space);
        }
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        
        meternmc = new MeterNormalModeCoordinate(coordinateDefinition, nm.getWaveVectorFactory().getWaveVectors());
        meternmc.setEigenVectors(nm.getEigenvectors(box));
        meternmc.setOmegaSquared(nm.getOmegaSquared(box));
        
        int coordNum = nm.getWaveVectorFactory().getWaveVectors().length*coordinateDim*2;
        hists = new AccumulatorHistogram[coordNum];
        DataSplitter splitter = new DataSplitter();
        DataPump pump = new DataPump(meternmc, splitter);

        
        for(int i = 0; i < coordNum; i++){
            hists[i] = new AccumulatorHistogram();
            splitter.setDataSink(i, hists[i]);
        }
        
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(1000);
        integrator.getEventManager().addListener(pumpListener);
        
        
        
        
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 1;
        int nA = 32;
        double density = 1.3;
        if (D == 1) {
            nA = 32;
            density = 0.5;
        }

        int nSteps = 10000;

        // parse arguments
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
//            simTime = Double.parseDouble(args[2]);
        }
        if (args.length > 3) {
            nA = Integer.parseInt(args[3]);
        }
        String filename = "normal_modes" + D + "D_"+nA+"_"+((int)(density*100))+"_cubic";
        if (args.length > 0) {
            filename = args[0];
        }
        
        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(nSteps + " time units");
        System.out.println("output data to " + filename);

        // construct simulation
        SimDegreeFreedom sim = new SimDegreeFreedom(Space.getInstance(D), nA, density);
        
        // start simulation
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.getController().actionPerformed();
        System.out.println("equilibration finished");
        sim.getController().reset();
        
        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();
        
        
        
        /* loops that
         *  -changes filename
         *  -changes accumulator histogram
         *  -changes histogram
         *  -calls actionPerformed
         *  
         */
//        int accumulatorLength = ;
//        for(int i = 0; i < ; i++){
//            
//            String filename = new String("hist_" + i);
//            wh = new WriteHistogram(filename);
//            wh.setHistogram();
//            wh.actionPerformed();
//            
//        }
        
        
        
        
        
        System.out.println("Fini.");
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public Boundary bdry;
    public Primitive primitive;
    public CoordinateDefinition coordinateDefinition;
    MeterNormalModeCoordinate meternmc;
    DoubleRange histRange;
    WaveVectorFactory waveVectorFactory;
    NormalModes nm;
    MCMoveAtomCoupled mcmove;
    AccumulatorHistogram[] hists;

}