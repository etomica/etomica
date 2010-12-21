package etomica.normalmode;

import etomica.action.BoxInflateAnisotropic;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.atom.DiameterHashByType;
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
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcpBaseCentered;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * NPT simulation for hard sphere solid using an MCMove that does coordinate
 * scaling with volume changes.
 */
public class HSDimerNPT extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public final PotentialMasterList potentialMaster;
    public final IntegratorMC integrator;
    public final SpeciesHSDimer species;
    public final IBox box;
    public final ActivityIntegrate activityIntegrate;
    public final CoordinateDefinitionHSDimer coordinateDefinition;

    public HSDimerNPT(Space space, int numMolecules, double rho, int[] nC) {
        super(space);
        potentialMaster = new PotentialMasterList(this, space);
        
        // Just initial setting
        double contB = 1.2;
        double contC = 1.4;
        double theta = Math.PI/2 - Math.asin(0.6/Math.sqrt(3.0));
        double a = Math.pow( (2.0/(contB*contC*rho))/Math.sin(theta), 1.0/3.0);        
        double b = contB*a;
        double c = contC*a;
        a = 1.2;
        b = Math.sqrt(3)*a;
        c = Math.sqrt(3)*a;
        
        System.out.println("a: " + a);
        System.out.println("b: " + b);
        System.out.println("c: " + c);
        
        double sigma = 1.0;
        
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(2.5);
        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        
        species = new SpeciesHSDimer(space, true);
        addSpecies(species);

        P2HardSphere p2 = new P2HardSphere(space, sigma, false);
        potentialMaster.addPotential(p2, new IAtomType[]{species.getDimerAtomType(), species.getDimerAtomType()});
    
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numMolecules);
        
        IVectorMutable[] boxDim = new IVectorMutable[3];
        boxDim[0] = space.makeVector(new double[]{nC[0]*a, 0.0, 0.0});
        boxDim[1] = space.makeVector(new double[]{0.0, nC[1]*b, 0.0});
        boxDim[2] = space.makeVector(new double[]{nC[2]*c*Math.cos(Math.PI/2+Math.asin(1.0/Math.sqrt(3.0))), 0.0, nC[2]*c*Math.sin(Math.PI/2+Math.asin(1.0/Math.sqrt(3.0)))});
        
		Primitive primitive = new PrimitiveMonoclinic(space, nC[0]*a, nC[1]*b, nC[2]*c, (Math.PI/2+Math.asin(1.0/Math.sqrt(3.0))));
		
		Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
        Basis basisHCPBase = new BasisHcpBaseCentered();
        Basis basis = new BasisBigCell(space, basisHCPBase, new int[]{nC[0], nC[1], nC[2]});
		box.setBoundary(boundary);
		
		for (int i=0; i<3; i++){
			System.out.println("boxDim["+i+"]: " + boxDim[i]);
		}
		
		System.out.println("0: " + box.getBoundary().getBoxSize().getX(0));
		System.out.println("1: " + box.getBoundary().getBoxSize().getX(1));
		System.out.println("2: " + box.getBoundary().getBoxSize().getX(2));
//		System.exit(1);
        coordinateDefinition = new CoordinateDefinitionHSDimer(this, box, primitive, basis, space);
        coordinateDefinition.setOrientationVectorCP2();
        coordinateDefinition.initializeCoordinates(new int[]{1,1,1});
        integrator.setBox(box);

        potentialMaster.getNeighborManager(box).reset();
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
       
        MCMoveMoleculeCoupled mcMove = new MCMoveMoleculeCoupled(potentialMaster,getRandom(),space);
        mcMove.setBox(box);
        mcMove.setDoExcludeNonNeighbors(true);
        //mcMove.setStepSize(0.01);
        integrator.getMoveManager().addMCMove(mcMove);

        MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
        rotate.setBox(box);
        //integrator.getMoveManager().addMCMove(rotate);
       
        
//        BoxInflateAnisotropic boxInfateAnis = new BoxInflateAnisotropic(box, space);
//        IVectorMutable scaleVec = space.makeVector(new double[]{0.9, 0.9, 0.9});
//        boxInfateAnis.setVectorScale(scaleVec);
//        boxInfateAnis.actionPerformed();
//        
//        System.exit(1);
        
        
        
        
        // using Carnahan-Starling EOS.  the pressure will be too high because
        // we have a solid, but OK.
        double eta = rho * Math.PI / 6;
        double den = 1-eta;
        double z = (1 + eta + eta*eta - eta*eta*eta) / (den*den*den);
        double p = z * rho;

        // hard coded pressure for rho=1.2
        p =10000;//23.3593;
        
        MCMove mcMoveVolume;
        if (false) {
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
        
        MCMoveVolumeMonoclinic mcMoveVolMonoclinic = new MCMoveVolumeMonoclinic(potentialMaster, getRandom(), space, p);
        mcMoveVolMonoclinic.setBox(box);
        mcMoveVolMonoclinic.setStepSize(0.001);
        ((MCMoveStepTracker)mcMoveVolMonoclinic.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(mcMoveVolMonoclinic);
        
        MCMoveVolumeMonoclinicAngle mcMoveVolMonoclinicAngle = new MCMoveVolumeMonoclinicAngle(potentialMaster, getRandom(), space, p, box);
        mcMoveVolMonoclinicAngle.setBox(box);
        mcMoveVolMonoclinicAngle.setStepSize(0.001);
        ((MCMoveStepTracker)mcMoveVolMonoclinicAngle.getTracker()).setNoisyAdjustment(true);
        integrator.getMoveManager().addMCMove(mcMoveVolMonoclinicAngle);
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HSMD3DParameters params = new HSMD3DParameters();
        ParseArgs parseArgs = new ParseArgs(params);
        parseArgs.parseArgs(args);
        
        int[] nC = params.nC;
        int numMolecules = nC[0]*nC[1]*nC[2]*2;
        HSDimerNPT sim = new HSDimerNPT(Space3D.getInstance(), numMolecules, params.rho, params.nC);
        
//        if(true){
//			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
//		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
//						
//			DiameterHashByType diameter = new DiameterHashByType(sim);
//			diameter.setDiameter(sim.species.getDimerAtomType(), 1.0);
//			simGraphic.getDisplayBox(sim.box).setDiameterHash(diameter);
//			
//			simGraphic.makeAndDisplayFrame("HS Dumbbell Crystal Structure");
//
//		    return;
//		    
//        }
        
        MeterVolume meterVolume = new MeterVolume();
        meterVolume.setBox(sim.box);
        DataFork volumeFork = new DataFork();
        DataPumpListener volumePump = new DataPumpListener(meterVolume, volumeFork, numMolecules);
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
			DiameterHashByType diameter = new DiameterHashByType(sim);
			diameter.setDiameter(sim.species.getDimerAtomType(), 1.0);
			graphic.getDisplayBox(sim.box).setDiameterHash(diameter);
			
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
        System.out.println("avg density "+numMolecules/vavg+" "+numMolecules/(vavg*vavg)*verr);
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
        public double rho = 1.2;
        public int[] nC = new int[]{4,4,4};
        public long numSteps = 10000000;
    }
}
