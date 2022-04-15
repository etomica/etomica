package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.data.*;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveGeneral;
import etomica.nbr.CriterionSimple;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialCuLREP;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.ElectronVolt;
import etomica.units.Kelvin;
import etomica.units.Pascal;
import etomica.units.Second;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

public class SimSolidPropsCuLREP extends Simulation {
    public PotentialMasterList potentialMaster;
    public IntegratorVelocityVerlet integrator;
    public SpeciesGeneral species;
    public Box box;
    public PotentialCuLREP potential;
    public DisplayBox display;
    public DisplayPlot plot;
    public MeterEnergy energy;
    public IDataInfo info2;
    public CoordinateDefinition coordinateDefinition;
    public Primitive primitive;
    public Boundary boundary;

    public SimSolidPropsCuLREP(int numAtoms, double density, double temperature, double dt, int mdHybSteps) {
        super(Space3D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.element(Copper.INSTANCE), true);
        addSpecies(species);

        potentialMaster = new PotentialMasterList(this, space);

        double L = Math.pow(4.0/density, 1.0/3.0);
        Vector[] cellDim = new Vector[3];
        cellDim[0] = Vector.of(new double[]{L, 0, 0});
        cellDim[1] = Vector.of(new double[]{0, L, 0});
        cellDim[2] = Vector.of(new double[]{0, 0, L});
        primitive = new PrimitiveGeneral(space, cellDim);

        int nC = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        int[] nCells = new int[]{nC,nC,nC};
        boundary = new BoundaryRectangularPeriodic(space,L*nC);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);
        Basis basis = new BasisCubicFcc();
        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        /**Hybrid MC*/
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);

        integrator.setThermostat(ThermostatType.HYBRID_MC);
        getController().addActivity(new ActivityIntegrate(integrator));
        integrator.setTimeStep(dt);
        integrator.setTemperature(temperature);
        integrator.setThermostatInterval(mdHybSteps);
        integrator.setIsothermal(true);
        integrator.setThermostatNoDrift(true);

        potential = new PotentialCuLREP(space);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[] {sphereType, sphereType });

        int cellRange = 4;
        double rcFacPotentialmaster = 1.0;
        double potentialMasterRc = potential.getRange()*rcFacPotentialmaster;
        potentialMaster.setRange(potentialMasterRc);
        potentialMaster.setCellRange(cellRange);
        potentialMaster.getNeighborManager(box).reset();
        potential.setRange(0.6*boundary.getBoxSize().getX(0), 0.6*boundary.getBoxSize().getX(0));
    }

    public static void main(String[] args) {
        Params parameters=new Params();
        ParseArgs.doParseArgs(parameters, args);
        int numAtoms = parameters.numAtoms;
        double density = parameters.density;
        double dt = parameters.dt;
        int mdHybSteps = parameters.mdHybSteps;
        int intervalData = parameters.intervalData;
        double temperature = parameters.temperature;
        long numSteps=parameters.numSteps;

        System.out.println(" "+numAtoms+" atoms , density: "+density+" atoms/A^3 , temperature: "+temperature + " K");
    	double a0 = Math.pow(4/density, 1.0/3.0);
    	System.out.println(" Volume/N: " + (1/density) + " A^3/atom ,  lattice constant: " + a0+ "  , Box size: " + ((int)Math.round(Math.pow(numAtoms/4, 1.0/3.0))*a0));
     
    	temperature = Kelvin.UNIT.toSim(temperature);

        final SimSolidPropsCuLREP sim = new SimSolidPropsCuLREP(numAtoms, density, temperature, dt, mdHybSteps);

        System.out.println("");
        System.out.println(" dt: " + Second.UNIT.fromSim(dt)*1E15 + " (fs)  " +"  ,  mdHybSteps: "+ mdHybSteps + " , intervalData: " + intervalData);

        int numBlocks = 100;
        int intervalMD = intervalData*mdHybSteps;
        long blockSize = numSteps/(numBlocks*intervalMD);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" "+ numBlocks + " blocks ; block size "+blockSize+" interval "+intervalMD);

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setBox(sim.box);
        double ULat = meterPE.getDataAsScalar();
        System.out.println(" uLat: " + ElectronVolt.UNIT.fromSim(ULat/numAtoms) + " (eV/atom)");
        
//        MeterPressure  meterP = new MeterPressure(sim.space);
//        meterP.setBox(sim.box);
//        meterP.setTemperature(0);//ZERO means NO ideal gas component!
//        meterP.setPotentialMaster(sim.potentialMaster);
//        double pLat = meterP.getDataAsScalar();
//        System.out.println(" pLat: " + (Pascal.UNIT.fromSim(pLat)/1.0e9) + " (GPa)");
//    
        
        MeterSolidProps meter = new MeterSolidProps(sim.getSpace(),meterPE,sim.potentialMaster,sim.coordinateDefinition, temperature, 0, 0, ULat, 0, false);
        
//        System.out.flush();
        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1000);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (allColors == null) {
                        allColors = new Color[768];
                        for (int i = 0; i < 256; i++) {
                            allColors[i] = new Color(255 - i, i, 0);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 256] = new Color(0, 255 - i, i);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 512] = new Color(i, 0, 255 - i);
                        }
                    }
                    return allColors[(2 * a.getLeafIndex()) % 768];
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame(" LREP ");

            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 10);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());
            return;
        }


/**Equilibaration*/
        final long startTime = System.currentTimeMillis();
        long Ninit = numSteps/5;
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Ninit));
        System.out.println(" "+ numSteps+" MD steps (or, " +(numSteps/mdHybSteps)+ " MC tials) after  " +Ninit+" MD steps of equil.");
        System.out.flush();

//        if(!false){
//            PotentialCalculationForceSum pcForce = new PotentialCalculationForceSum();
//            AtomLeafAgentManager.AgentSource<IntegratorVelocityVerlet.MyAgent> atomAgentSource = new AtomLeafAgentManager.AgentSource<IntegratorVelocityVerlet.MyAgent>() {
//                public MyAgent makeAgent(IAtom a, IBox agentBox) {
//                	IntegratorVelocityVerlet.MyAgent agent = new IntegratorVelocityVerlet.MyAgent(sim.space);
//                    return agent;
//                }
//                @Override
//                public void releaseAgent(MyAgent agent, IAtom atom, IBox agentBox) {
//                }
//            };            
//            AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent> atomAgentManager = new AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent>(atomAgentSource , sim.box , IntegratorVelocityVerlet.MyAgent.class);
//            pcForce.setAgentManager(atomAgentManager);
//            IteratorDirective id = new IteratorDirective();
//            id.includeLrc = false;
//
//            IVectorMutable F  = sim.space.makeVector();
//
//            
//            MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(sim.potentialMaster);
//            meterPE2.setBox(sim.box);
//            int N = 200;
//            IVectorMutable[] dr = new IVectorMutable[numAtoms];
//            IAtom atom;
//	        for(int i=0;i<numAtoms;i++){
//	        	dr[i] = sim.space.makeVector();
//	        	atom = sim.box.getLeafList().getAtom(i);
//	        	IVector Ri = sim.coordinateDefinition.getLatticePosition(sim.box.getLeafList().getAtom(i));
//	        	IVector ri = atom.getPosition();
//               	dr[i].Ev1Mv2(ri, Ri);
//	        }
//
//            for(int n=0;n<=N;n++){
//	        	double frac = n/((double)N);
//    	        for(int i=0;i<numAtoms;i++){
//    	        	atom = sim.box.getLeafList().getAtom(i);
//    	        	IVector Ri = sim.coordinateDefinition.getLatticePosition(atom);
//    	        	IVectorMutable ri = atom.getPosition();
//    	        	ri.E(Ri);
//    	        	ri.PEa1Tv1(frac, dr[i]);
//    	        }
//                sim.potentialMaster.calculate(sim.box, id, pcForce);
//    	        double Fdr = 0;
//    	        for(int i=0;i<numAtoms;i++){
//    	        	atom = sim.box.getLeafList().getAtom(i);
//    	            F.E(((IntegratorVelocityVerlet.MyAgent)atomAgentManager.getAgent(atom)).force);
//    	            Fdr += frac*F.dot(dr[i]); 
//    	        }
//    	         pcForce.reset(); 
//
//    	        double dU = meterPE2.getDataAsScalar() - ULat;
//    	        
//    	        System.out.println(frac + "   "+ dU + "   " +0.5*Fdr +"   "+ (dU+0.5*Fdr));
//    	        
//    	        
//    	        for(int i=0;i<numAtoms;i++){
//    	        	IVector Ri = sim.coordinateDefinition.getLatticePosition(sim.box.getLeafList().getAtom(i));
//    	        	IVectorMutable ri = sim.box.getLeafList().getAtom(i).getPosition();
//    	        	ri.Ev1Pv2(Ri, dr[i]);
//    	        }
//
//           }
//
//            
//            
//      } // END debug
//
        
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPump = new DataPumpListener(meter, accumulator, intervalMD);
        sim.integrator.getEventManager().addListener(accumulatorPump);

        
/***********/

//RUN...        
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));






//        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(sim.potentialMaster);
//        meterPE2.setBox(sim.box);
//        int N = 200;
//        IVectorMutable[] dri = new IVectorMutable[numAtoms];
//        for(int n=1;n<=N;n++){
//	        for(int i=0;i<numAtoms;i++){
//	        	dri[i] = sim.space.makeVector();
//	        	IVector Ri = sim.coordinateDefinition.getLatticePosition(sim.box.getLeafList().getAtom(i));
//	        	IVectorMutable ri = sim.box.getLeafList().getAtom(i).getPosition();
//               	dri[i].Ev1Mv2(ri, Ri);
//	        	double frac = n/((double)N);
//	        	ri.E(Ri);
//	        	ri.PEa1Tv1(frac, dri[i]);
//	        }
//	        double dU = meterPE2.getDataAsScalar();
//	        
//	        System.out.println(n + " "+ (dU-ULat)/numAtoms);
//	        
//	        
//	        for(int i=0;i<numAtoms;i++){
//	        	IVector Ri = sim.coordinateDefinition.getLatticePosition(sim.box.getLeafList().getAtom(i));
//	        	IVectorMutable ri = sim.box.getLeafList().getAtom(i).getPosition();
//	        	ri.Ev1Pv2(Ri, dri[i]);
//	        }
//
//       }
        
        
        
        
        
        
        double accRatio = sim.integrator.getHybridAcceptance();
        System.out.println("");
        System.out.println(" accRatio: " + accRatio);

        DataGroup data = (DataGroup)accumulator.getData();
        IData dataAvg  = data.getData(accumulator.AVERAGE.index);
        IData dataErr  = data.getData(accumulator.ERROR.index);
        IData dataCorr = data.getData(accumulator.BLOCK_CORRELATION.index);

        double dU   = dataAvg.getValue(0);
//        double Ur   = dataAvg.getValue(1); 
        System.out.println();
        System.out.println("*******************************************************************************************************");
        System.out.println(" conv.  uAnhC(eV/atom): " + ElectronVolt.UNIT.fromSim(dU-1.5*(numAtoms-1)*temperature)/numAtoms + "  +/- " + ElectronVolt.UNIT.fromSim(dataErr.getValue(0))/numAtoms + "  corr: "+ dataCorr.getValue(0));
//        System.out.println(" mapped uAnhM(eV/atom): " +ElectronVolt.UNIT.fromSim(Ur)/numAtoms + "  +/- " 
//                + ElectronVolt.UNIT.fromSim(dataErr.getValue(1))/numAtoms + "  corr: "+ dataCorr.getValue(1));
//        System.out.println("*******************************************************************************************************\n");


        
//        System.out.println("difference(eV): " + ElectronVolt.UNIT.fromSim((ULat+dU)/numAtoms-(ULat+ 1.5*(numAtoms-1)*temperature + Ur)/numAtoms));
//        System.out.println("       ratiodU: " + dataErr.getValue(0)/dataErr.getValue(1));
        System.out.println();
        long endTime = System.currentTimeMillis();
//        System.out.println("time(hrs): " + (endTime - startTime)/1000.0/3600.0);
        System.out.println("time(sec): " + (endTime - startTime)/1000.0);

    }    
    public static class Params extends ParameterBase{
    	public int numAtoms = 500;
        public double density = 0.08502338387498792;
//        public double density = 0.1062792298437349;
        public double temperature = 1086.2160000000001;
        public double dt =  0.004;
        public long numSteps = 100000;
        public int  mdHybSteps  = 5;
        public int intervalData = 1;
     }
}