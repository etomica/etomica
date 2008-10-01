package etomica.models.oneDHardRods;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.data.DataPump;
import etomica.data.DataSource;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.NormalModes1DHR;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;
import etomica.virial.overlap.IntegratorOverlap;

public class SimOverlapAB extends Simulation {
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "SimOverlapAB";
    CoordinateDefinition coordinateDefinition;
    Primitive primitive;
    Basis basis;
    int[] nCells;
    NormalModes1DHR nm;
    int affectedWV;
    double alpha;       //adjustable parameter
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    ActivityIntegrate activityIntegrate;
    
    
    IntegratorMC[] integrators;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] pumps;
    public DataSource[] meters;
    public IBox boxTarget, boxRef;
    public Boundary boundaryTarget, boundaryRef;

    
    public SimOverlapAB(Space _space, int numAtoms, double density, double 
            temperature, String filename, double harmonicFudge){
        super(_space, true);
        
        //Set up some of the joint stuff
        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        
        integrators = new IntegratorMC[2];
        pumps = new DataPump[2];
        meters = new DataSource[2];
        accumulators = new AccumulatorVirialOverlapSingleAverage[2];
        
        
        
        //Set up target system
        
        
        
        
        
        
        
        
        
        //Set up reference system
        
        
        
        
        
        
        //Set up the rest of the joint stuff
    }
    
    
    
    public static void main(String args[]){
        
        
    }
    
    public static class SimOverlapABParam extends ParameterBase {
        public int numAtoms = 32;
        public double density = 0.5;
        public int D = 1;
        public long numSteps = 1000000;
        public double harmonicFudge = 1.0;
        public String filename = "HR1D_";
        public double temperature = 1.0;
        public int affectedWV = 16;
    }
    
}
