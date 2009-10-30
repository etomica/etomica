package etomica.models.oneDHardRods;

import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.IPotentialMaster;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.units.Null;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * Uses a different box than the main simulation, to assume an extra mode & rod 
 * is added/
 * 
 * @author cribbin
 *
 */
public class MeterDifferentImage extends DataSourceScalar {

    public int nInsert, counter;
    private MeterPotentialEnergy meterPE;
    private CoordinateDefinition coordinateDefinition, oldCD;
    private int coordinateDim;
    private double eigenVectors[][][];
    private IVectorMutable[] waveVectors;
    private double[] realT, imagT;
    private double[][] uOld, omegaSquared;
    protected double temperature;
    private double[] uNow, deltaU;
    private double[] waveVectorCoefficients;
    
    private IBox box;
    private int numAtoms;
    private Boundary bdry;
    
    
    public MeterDifferentImage(String string, IPotentialMaster potentialMaster, 
            IPotential potential, int numSimAtoms, double density, Simulation sim,
            Primitive simPrimitive, Basis simBasis, CoordinateDefinition simCD){
        super(string, Null.DIMENSION);
        
        oldCD = simCD;
        
        numAtoms = numSimAtoms + 1;
        coordinateDim = numAtoms * 2;
        box = new Box(sim.getSpace());
        box.setNMolecules(sim.getSpeciesManager().getSpecies(0), numAtoms); 
        sim.addBox(box);
        
        potential.setBox(box);
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        
        bdry = new BoundaryRectangularPeriodic(sim.getSpace(), numAtoms/density);
        box.setBoundary(bdry);
        
        int[] nCells = new int[]{numAtoms};
        coordinateDefinition = new CoordinateDefinitionLeaf(box, simPrimitive, 
                simBasis, sim.getSpace());
        coordinateDefinition.initializeCoordinates(nCells);
        
        realT = new double[coordinateDim];
        imagT = new double[coordinateDim];
        deltaU = new double[coordinateDim];
        

    }
    
    public double getDataAsScalar() {

        //Calculate normal mode coordinates of simulation system.
        
        
        
        //Assign the last normal mode coordinate from the Gaussian distribution
        
        
        
        //Assign that set of normal mode coordinates to the meter's system.
        
        
        
        //Calculate the positions for the meter's system
        
        
        
        //Check for overlap
        
        
        
        
        
        double energy = meterPE.getDataAsScalar();
        if(Double.isInfinite(energy)) {
            return 0;
        } else {
            return 1;
        }
    }

    
    
}
