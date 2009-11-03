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
    private CoordinateDefinition cDef, simCDef;
    private int cDim, simCDim;
    private double eigenVectors[][][];
    private IVectorMutable[] waveVectors, simWaveVectors;
    private double[] realT, imagT, simRealT, simImagT;
    private double[][] uOld, omegaSquared;
    protected double temperature;
    private double[] uNow, deltaU;
    private double[] waveVectorCoefficients;
    
    private IBox box;
    private int numAtoms;
    private Boundary bdry;
    
    
    public MeterDifferentImage(String string, IPotentialMaster potentialMaster, 
            IPotential potential, int numSimAtoms, double density, Simulation sim,
            Primitive simPrimitive, Basis simBasis, CoordinateDefinition simCD,
            IVectorMutable[] simWV){
        super(string, Null.DIMENSION);
        
        simWaveVectors = simWV;
        this.simCDef = simCD;
        simCDim = simCD.getCoordinateDim();
        
        numAtoms = numSimAtoms + 1;
        box = new Box(sim.getSpace());
        box.setNMolecules(sim.getSpeciesManager().getSpecies(0), numAtoms); 
        sim.addBox(box);
        
        potential.setBox(box);
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        
        bdry = new BoundaryRectangularPeriodic(sim.getSpace(), numAtoms/density);
        box.setBoundary(bdry);
        
        int[] nCells = new int[]{numAtoms};
        cDef = new CoordinateDefinitionLeaf(box, simPrimitive, 
                simBasis, sim.getSpace());
        cDef.initializeCoordinates(nCells);
        cDim = cDef.getCoordinateDim();
        
        realT = new double[cDim];
        imagT = new double[cDim];
        simRealT = new double[simCDim];
        simImagT = new double[simCDim];

    }
    
    public double getDataAsScalar() {

        //Store old positions
        BasisCell[] simCells = simCDef.getBasisCells();
        BasisCell cell = simCells[0];
        uOld = new double[simCells.length][simCDim];
        for(int iCell = 0; iCell < simCells.length; iCell++){
            //store old positions.
            uNow = simCDef.calcU(simCells[iCell].molecules);
            System.arraycopy(uNow, 0, uOld[iCell], 0, cDim);
        }
        
        //Calculate normal mode coordinates of simulation system.
        double[] simRealCoord = new double[cDim - 1];
        double[] simImagCoord = new double[cDim - 1];
        
        for (int wvcount = 0; wvcount < waveVectors.length; wvcount++){
            simCDef.calcT(waveVectors[wvcount], simRealT, simImagT);
            for (int iCell = 0; iCell < simCells.length; iCell++){
                cell = simCells[iCell];
                
                //Calculate the contributions to the current position of this mode
                double kr = simWaveVectors[wvcount].dot(cell.cellPosition);
            }
        }
        
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
