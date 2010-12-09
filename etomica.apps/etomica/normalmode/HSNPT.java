package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterVolume;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * NPT simulation for hard sphere solid using an MCMove that does coordinate
 * scaling with volume changes.
 */
public class HSNPT extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public final PotentialMasterList potentialMaster;
    public final IntegratorMC integrator;
    public final SpeciesSpheresMono species;
    public final IBox box;
    public final ActivityIntegrate activityIntegrate;
    public final CoordinateDefinition coordinateDefinition;

    public HSNPT(Space _space, int numAtoms, double rho) {
        super(_space);
        potentialMaster = new PotentialMasterList(this, space);
        
        double neighborRangeFac = 1.4;
        double sigma = 1.0;
        double l = Math.pow(numAtoms / rho, 1.0/3.0);
        potentialMaster.setCellRange(1);
        potentialMaster.setRange(neighborRangeFac*sigma);
        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        IAtomType type1 = species.getLeafType();

        P2HardSphere p2 = new P2HardSphere(space, sigma, false);
        potentialMaster.addPotential(p2, new IAtomType[]{type1, type1});

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{l,l,l}));
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        coordinateDefinition = new CoordinateDefinitionLeaf(box, new PrimitiveCubic(space, l/n), new BasisCubicFcc(), space);
        coordinateDefinition.initializeCoordinates(new int[]{n,n,n});
        integrator.setBox(box);

        potentialMaster.getNeighborManager(box).reset();
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        MCMoveAtomCoupled mcMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        mcMove.setPotential(p2);
        integrator.getMoveManager().addMCMove(mcMove);

        // using Carnahan-Starling EOS.  the pressure will be too high because
        // we have a solid, but OK.
        double eta = rho * Math.PI / 6;
        double den = 1-eta;
        double z = (1 + eta + eta*eta - eta*eta*eta) / (den*den*den);
        double p = z * rho;

        // hard coded pressure for rho=1.2
        p = 23.3593;

        MCMove mcMoveVolume;
        if (true) {
            // fancy move
            mcMoveVolume = new MCMoveVolumeSolid(potentialMaster, coordinateDefinition, getRandom(), space, p);
            ((MCMoveVolumeSolid)mcMoveVolume).setTemperature(1.0);
        }
        else {
            // standard move
            mcMoveVolume = new MCMoveVolume(potentialMaster, getRandom(), space, p);
        }
        ((MCMoveStepTracker)mcMoveVolume.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(mcMoveVolume);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD3DParameters params = new HSMD3DParameters();
        ParseArgs parseArgs = new ParseArgs(params);
        parseArgs.parseArgs(args);
        HSNPT sim = new HSNPT(Space3D.getInstance(), params.numAtoms, params.rho);
        
        MeterVolume meterVolume = new MeterVolume();
        meterVolume.setBox(sim.box);
        DataFork volumeFork = new DataFork();
        DataPumpListener volumePump = new DataPumpListener(meterVolume, volumeFork, params.numAtoms);
        sim.integrator.getEventManager().addListener(volumePump);
        AccumulatorAverageCollapsing volumeAvg = new AccumulatorAverageCollapsing();
        volumeFork.addDataSink(volumeAvg);
        
//        MeterDisplacement meterDisplacement = new MeterDisplacement(sim.space, sim.coordinateDefinition, 0.001);
//        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterDisplacement, params.numAtoms));

//        MeterDisplacementRMS meterDisplacementRMS = new MeterDisplacementRMS(sim.space, sim.coordinateDefinition, 0.001);
//        AccumulatorAverageFixed displacementAvg = new AccumulatorAverageFixed();
//        DataPumpListener displacementAvgPump = new DataPumpListener(meterDisplacementRMS, displacementAvg, params.numAtoms);
//        sim.integrator.getEventManager().addListener(displacementAvgPump);
//
//        MeterDisplacementMax meterDisplacementMax = new MeterDisplacementMax(sim.space, sim.coordinateDefinition, 0.001);
//        AccumulatorAverageFixed displacementMax = new AccumulatorAverageFixed();
//        DataPumpListener displacementMaxPump = new DataPumpListener(meterDisplacementMax, displacementMax, params.numAtoms);
//        sim.integrator.getEventManager().addListener(displacementMaxPump);
//
//        MeterMaxExpansion meterMaxExpansion = new MeterMaxExpansion(sim.space, sim.box, sim.potentialMaster.getNeighborManager(sim.box));
//        AccumulatorAverageFixed maxExpansionAvg = new AccumulatorAverageFixed();
//        DataPumpListener maxExpansionPump = new DataPumpListener(meterMaxExpansion, maxExpansionAvg, params.numAtoms);
//        sim.integrator.getEventManager().addListener(maxExpansionPump);

        if (true) {
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, sim.getSpace(), sim.getController());

            AccumulatorHistory densityHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//            densityHistory.setPushInterval();
            volumeFork.addDataSink(densityHistory);
            DisplayPlot densityPlot = new DisplayPlot();
            densityHistory.setDataSink(densityPlot.getDataSet().makeDataSink());
            densityPlot.setLabel("volume");
            graphic.add(densityPlot);
            
            graphic.makeAndDisplayFrame();
            
            return;
        }
        sim.activityIntegrate.setMaxSteps(params.numSteps/10);
        sim.activityIntegrate.actionPerformed();
        volumeAvg.reset();
        System.out.println("equilibration finished");
        sim.activityIntegrate.setMaxSteps(params.numSteps);
        sim.activityIntegrate.actionPerformed();
//        try {
//            FileWriter fw = new FileWriter("disp"+params.rho+".dat");
//            IData data = meterDisplacement.getData();
//            IData xData = meterDisplacement.getIndependentData(0);
//            for (int i=0; i<data.getLength(); i++) {
//                fw.write(xData.getValue(i)+" "+data.getValue(i)+"\n");
//            }
//            fw.close();
//        }
//        catch (IOException e) {
//            throw new RuntimeException(e);
//        }
        
        double vavg = volumeAvg.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
        double verr = volumeAvg.getData().getValue(AccumulatorAverage.StatType.ERROR.index);
        double vstdev = volumeAvg.getData().getValue(AccumulatorAverage.StatType.STANDARD_DEVIATION.index);
        System.out.println("avg density "+params.numAtoms/vavg+" "+params.numAtoms/(vavg*vavg)*verr);
        System.out.println("avg volume "+vavg+" stdev "+vstdev);

//        double davg = displacementAvg.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
//        double dstdev = displacementAvg.getData().getValue(AccumulatorAverage.StatType.STANDARD_DEVIATION.index);
//        System.out.println("displacement avg "+davg+" stdev "+dstdev);
//
//        double dmaxavg = displacementMax.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
//        double dmaxstdev = displacementMax.getData().getValue(AccumulatorAverage.StatType.STANDARD_DEVIATION.index);
//        System.out.println("displacement max avg "+dmaxavg+" stdev "+dmaxstdev);
//
//        double emaxavg = maxExpansionAvg.getData().getValue(AccumulatorAverage.StatType.AVERAGE.index);
//        double emaxstdev = maxExpansionAvg.getData().getValue(AccumulatorAverage.StatType.STANDARD_DEVIATION.index);
//        System.out.println("max expansion avg "+emaxavg+" stdev "+emaxstdev);
    }
    
    public static class HSMD3DParameters extends ParameterBase {
        public int numAtoms = 500;
        public double rho = 1.25;
        public long numSteps = 10000000;
    }
}
