package etomica.modules.dcvgcmd;
import etomica.ConfigurationLattice;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceGroup;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfile;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2WCA;
import etomica.space.BoundaryRectangularSlit;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.Kelvin;

public class DCVGCMD extends Simulation {

	public IntegratorDCVGCMD integratorDCV;
	public P2WCA potential;
	public P2WCA potential1;
	public P1LJWCAWall potentialwall;
	public P1LJWCAWall potentialwall1;
	public SpeciesSpheresMono species;
	public SpeciesSpheresMono species1;
	public Phase phase;
	public DataSourceGroup fluxMeters;
	public MeterTemperature thermometer;
	public MeterNMolecules density1;
	public MeterNMolecules density2;
	public MeterProfile profile1;
	public MeterProfile profile2;
	public AccumulatorAverage accumulator1;
	public AccumulatorAverage accumulator2;
    public AccumulatorAverage fluxAccumulator;
	
	
 //Constructor
   public DCVGCMD() {
     //Instantiate classes
      super(new Space3D());
      Default.ATOM_MASS = 40.;
      Default.ATOM_SIZE = 3.41;
      Default.POTENTIAL_WELL = 119.8;
      //Default.makeLJDefaults();

      //Default.BOX_SIZE = 14.0;

        potential = new P2WCA(space);
		P2WCA potential11 = new P2WCA(space);
		potential1 = new P2WCA(space);
		potentialwall = new P1LJWCAWall(space);
		potentialwall1 = new P1LJWCAWall(space);
		species = new SpeciesSpheresMono(this);
		species1 = new SpeciesSpheresMono(this);
		potentialMaster.setSpecies(potential, new Species[] {species, species});
		potentialMaster.setSpecies(potential1, new Species[] {species1, species});
		potentialMaster.setSpecies(potential11, new Species[] {species1, species1});
		potentialMaster.setSpecies(potentialwall, new Species[] {species});
		potentialMaster.setSpecies(potentialwall1, new Species[] {species1});
		species.setNMolecules(8);
		species1.setNMolecules(8);
		phase = new Phase(this);
		integratorDCV = new IntegratorDCVGCMD(potentialMaster, species, species1);
        Default.AUTO_REGISTER = false;
		final IntegratorVelocityVerlet integrator = new IntegratorVelocityVerlet(potentialMaster, space);
        final IntegratorMC integratorMC = new IntegratorMC(potentialMaster);
        Default.AUTO_REGISTER = true;
	    integratorDCV.addPhase(phase);
	    
	    ActivityIntegrate activityIntegrate = new ActivityIntegrate(integratorDCV);
	    getController().addAction(activityIntegrate);
	    
	    //make MC integrator next
	    integratorDCV.setIntegrators(integratorMC, integrator);
        integratorDCV.setTemperature(Kelvin.UNIT.toSim(500.));
	    integrator.setIsothermal(true);
            //integrator.setSleepPeriod(1);
        integrator.setTimeStep(0.05);
            //integrator.setInterval(10);
        activityIntegrate.setDoSleep(true);	
        phase.setBoundary(new BoundaryRectangularSlit(space, 2)); 
        phase.boundary().setDimensions(new Vector3D(40,40,80));
        // Crystal crystal = new Crystal(new PrimitiveTetragonal(space, 20, 40),new BasisMonatomic(3));
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc()); 
        phase.setConfiguration(config);

        MyMCMove[] moves = integratorDCV.mcMoves();
		MeterFlux meterFlux0 = new MeterFlux(moves[0], integratorDCV);
		MeterFlux meterFlux1 = new MeterFlux(moves[1], integratorDCV);
		MeterFlux meterFlux2 = new MeterFlux(moves[2], integratorDCV);
		MeterFlux meterFlux3 = new MeterFlux(moves[3], integratorDCV);
        meterFlux0.setPhase(phase);
        meterFlux1.setPhase(phase);
        meterFlux2.setPhase(phase);
        meterFlux3.setPhase(phase);
		fluxMeters = new DataSourceGroup(new MeterFlux[] {meterFlux0, meterFlux1, meterFlux2, meterFlux3});
		fluxAccumulator = new AccumulatorAverage();
        DataPump fluxPump = new DataPump(fluxMeters, fluxAccumulator);
        IntervalActionAdapter fluxInterval = new IntervalActionAdapter (fluxPump, integratorDCV);
        
        thermometer = new MeterTemperature();
        thermometer.setPhase(phase);
        
        density1 = new MeterNMolecules();
		density2 = new MeterNMolecules();
		density1.setPhase(phase);
    	density2.setPhase(phase);
        density1.setSpecies(species);
    	density2.setSpecies(species1);
    	profile1 = new MeterProfile(space);
		profile2 = new MeterProfile(space);
		profile1.setMeter(density1);
    	profile1.setPhase(phase);
    	profile1.setProfileVector(new Vector3D(0.0,0.0,1.0));
    	profile2.setMeter(density2);
    	profile2.setPhase(phase);
    	profile2.setProfileVector(new Vector3D(0.0,0.0,1.0));
    	
    	accumulator1 = new AccumulatorAverage();
        DataPump profile1pump = new DataPump(profile1, accumulator1);
        new IntervalActionAdapter (profile1pump, integratorDCV);
    	
    	accumulator2 = new AccumulatorAverage();
		DataPump profile2pump = new DataPump(profile2, accumulator2);
    	IntervalActionAdapter interval2 = new IntervalActionAdapter (profile2pump, integratorDCV);
    	
        integrator.addIntervalListener(new PhaseImposePbc(phase));
   } //End of constructor

} //End of DCVGCMD class
