package etomica.normalmode;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPumpListener;
import etomica.data.DataSplitter;
import etomica.data.DataTag;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.DeviceBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;
import etomica.util.HistogramCollapsing;

/**
 * Applet to illustrate HTTP method
 * 
 * @author Tai Boon Tan
 */
public class SimOverlapSoftSphereTPApplet extends Simulation {

    public SimOverlapSoftSphereTPApplet(Space _space, int numAtoms,  double density, double temperature, double otherTemperature, int exponent, double rc) {
        super(_space);
        
        potentialMaster = new PotentialMasterList(this, space);
        
        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        atomMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        
    
        double L = Math.pow(4.0/density, 1.0/3.0);
        double nbrDistance = L / Math.sqrt(2);
        int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
        primitive = new PrimitiveCubic(space, n*L);
        
        nCells = new int[]{n,n,n};
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);
    
        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1,1,1});

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
     	if(potentialMaster instanceof PotentialMasterList){
			potential = new P2SoftSphericalTruncated(space, potential, rc);
		
		} else {
			potential = new P2SoftSphericalTruncatedShifted(space, potential, rc);
			
		}
        atomMove.setPotential(potential);
        IAtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new IAtomType[] {sphereType, sphereType });

        /*
         *  1-body Potential to Constraint the atom from moving too far
         *  	away from its lattice-site
         *  
         */

        p1Constraint = new P1ConstraintNbr(space, nbrDistance, this);
        p1Constraint.initBox(box);
        atomMove.setConstraint(p1Constraint);

        potentialMaster.lrcMaster().setEnabled(false);
    
        integrator.setBox(box);

		if (potentialMaster instanceof PotentialMasterList) {
            int cellRange = 7;
            ((PotentialMasterList)potentialMaster).setRange(rc);
            ((PotentialMasterList)potentialMaster).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMaster).getNeighborManager(box).reset();
            int potentialCells = ((PotentialMasterList)potentialMaster).getNbrCellManager(box).getLattice().getSize()[0];
            if (potentialCells < cellRange*2+1) {
                throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
            }
		}
        
		latticeEnergy = meterPE.getDataAsScalar();

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        
        if (potentialMaster instanceof PotentialMasterList) {
            // extend potential range, so that atoms that move outside the truncation range will still interact
            // atoms that move in will not interact since they won't be neighbors
            ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundary.getBoxSize().getX(0));
        }
    }
    
    public void setT1(double temperature){
    	integrator.setTemperature(temperature);
    	meter.setTemperature(temperature);
    }
    
    public double getT1(){
    	return integrator.getTemperature();
    }
    
    public void setT2(double temperature){
    	meter.setOtherTemperature(temperature);
    }
    
    public double getT2(){
    	return meter.otherTemperature;
    }
    
//    public class Applet extends SimulationGraphic{
//    	public Applet(SimOverlapSoftSphereTPApplet sim, Space space){
//    		
//    	}
//    }
    
    public static void main(String[] args) {
        //set up simulation parameters
        double density = 1.1964;
        int exponentN = 12;
        final int numMolecules = 32;
        double temperature = 0.5;
        double rc = 2.2;
        double otherTemperature = 1.0;

    	SimOverlapSoftSphereTPApplet sim = new SimOverlapSoftSphereTPApplet(Space.getInstance(3), numMolecules, density, temperature, otherTemperature, exponentN, rc);
   
        sim.meter = new MeterBoltzmannHTTP(sim.potentialMaster, sim.species, sim.space, sim);
        sim.meter.setCoordinateDefinition(sim.coordinateDefinition);
        sim.meter.setLatticeEnergy(sim.latticeEnergy);
        sim.meter.setTemperature(temperature);
        sim.meter.setOtherTemperature(otherTemperature);
        sim.meter.setConstraint(sim.p1Constraint);
            
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, sim.space, sim.getController());
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(80));

        DataSplitter  histogramSplitter = new DataSplitter();
        DataPumpListener boltzmannPump = new DataPumpListener(sim.meter, histogramSplitter, 50);
        sim.integrator.getEventManager().addListener(boltzmannPump);
            
        DisplayPlot boltzmannPlot = new DisplayPlot();
        final AccumulatorHistogram[] accumulatorHistogram = new AccumulatorHistogram[5];
        
        for(int i=0; i<5; i++){
        	accumulatorHistogram[i] = new AccumulatorHistogram(new HistogramCollapsing());
        	histogramSplitter.setDataSink(i, accumulatorHistogram[i]);
        	accumulatorHistogram[i].addDataSink(boltzmannPlot.getDataSet().makeDataSink());
        	accumulatorHistogram[i].setPushInterval(100);
        }
        
        boltzmannPlot.setLabel("Reduced Energy Histogram");
        
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[0].getTag()}, "beta1*u");
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[1].getTag()}, "beta2*u");
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[2].getTag()}, "beta2*u'");
      
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[3].getTag()}, "(beta2 - beta1)*u");
        boltzmannPlot.setLegend(new DataTag[]{accumulatorHistogram[4].getTag()}, "(beta2*u' - beta1*u");

        simGraphic.add(boltzmannPlot);
        
        DeviceBox t1Box = new DeviceBox();
        t1Box.setModifier(new ModifierGeneral(sim, "t1"));
        t1Box.setLabel("Temperature 1");
        
        DeviceBox t2Box = new DeviceBox();
        t2Box.setModifier(new ModifierGeneral(sim, "t2"));
        t2Box.setLabel("Temperature 2");
        
        simGraphic.add(t1Box);
        simGraphic.add(t2Box);

        simGraphic.getController().getResetAveragesButton().setPostAction(new IAction(){
        	public void actionPerformed(){
        		for (int i=0; i<5; i++){
        			accumulatorHistogram[i].reset();
        		}
        	}
        });
        
        simGraphic.makeAndDisplayFrame("HTTP Method");
  
    }
    
    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public AccumulatorAverageFixed accumulator;
    public DataPumpListener accumulatorPump;
    protected MCMoveAtomCoupled atomMove;
    protected PotentialMaster potentialMaster;
    protected double latticeEnergy;
    protected SpeciesSpheresMono species;
    protected P1ConstraintNbr p1Constraint;
    protected CoordinateDefinitionLeaf coordinateDefinition;
    protected MeterBoltzmannHTTP meter;
    
}
