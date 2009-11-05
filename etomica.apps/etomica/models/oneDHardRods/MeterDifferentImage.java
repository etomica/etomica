package etomica.models.oneDHardRods;

import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
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
    private IVectorMutable[] waveVectors, simWaveVectors;
    private double[] rRand, iRand, realT, imagT, simRealT, simImagT;
    private double[][] uOld, stdDev;
    protected double temperature;
    private double[] uNow, deltaU;
    private double[] wvCoeff, simWVCoeff;
    private double[][][] eigenVectors, simEigenVectors;

    protected final IRandom random;
    private IBox box;
    private int numAtoms;
    private Boundary bdry;
    
    
    public MeterDifferentImage(String string, IPotentialMaster potentialMaster, 
            int numSimAtoms, double density, Simulation sim,
            Primitive simPrimitive, Basis simBasis, CoordinateDefinition simCD,
            IVectorMutable[] simWV, double[][][] simEV){
        super(string, Null.DIMENSION);
        
        simWaveVectors = simWV;
        this.simCDef = simCD;
        simCDim = simCD.getCoordinateDim();
        simEigenVectors = simEV;
        this.random = sim.getRandom();
        
        numAtoms = numSimAtoms + 1;
        box = new Box(sim.getSpace());
        box.setNMolecules(sim.getSpecies(0), numAtoms); 
        sim.addBox(box);
        
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

        BasisCell[] simCells = simCDef.getBasisCells();
        BasisCell[] cells = cDef.getBasisCells();
        BasisCell cell = simCells[0];
        
        //Calculate normal mode coordinates of simulation system.
        double[] realCoord = new double[cDim];
        double[] imagCoord = new double[cDim];
        
        for (int wvcount = 0; wvcount < simWVCoeff.length; wvcount++){
            simCDef.calcT(simWaveVectors[wvcount], simRealT, simImagT);
            for (int i = 0; i < simCDim; i++){
                realCoord[wvcount] = 0.0;
                imagCoord[wvcount] = 0.0;
                for (int j = 0; j < simCDim; j++){
                    realCoord[wvcount] += simEigenVectors[wvcount][i][j] * simRealT[j];
                    imagCoord[wvcount] += simEigenVectors[wvcount][i][j] * simImagT[j];
                }
            }
        }
        
        //Create the last normal mode coordinate from the Gaussian distribution
        double sqrtT = Math.sqrt(temperature);
        for (int j = 0; j < cDim; j++) {
            if (stdDev[waveVectors.length][j] == 0) {continue;}
            // generate real and imaginary parts of random normal-mode coordinate 
            double realGauss = random.nextGaussian() * sqrtT;
            double imagGauss = random.nextGaussian() * sqrtT;
            if (wvCoeff[waveVectors.length] == 0.5) {imagGauss = 0;}
            realCoord[waveVectors.length] = realGauss * stdDev[waveVectors.length][j];
            imagCoord[waveVectors.length] = imagGauss * stdDev[waveVectors.length][j];
        }
            
        //Calculate the positions for the meter's system
        for (int wvcount = 0; wvcount < waveVectors.length; wvcount++){
            for (int iCell = 0; iCell < cells.length; iCell++){
                cell = cells[iCell];
                
                //Calculate the change in positions.
                double kR = waveVectors[wvcount].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for (int i = 0; i < cDim; i++){
                    for (int j = 0; j < cDim; j++){
                        deltaU[j] += wvCoeff[wvcount] * eigenVectors[wvcount][i][j] 
                            * 2.0 * (realCoord[i] * coskR - imagCoord[i] * sinkR);
                    }
                }
            }
        }
        
        
        //Check for overlap
        double energy = meterPE.getDataAsScalar();
        
        if(Double.isInfinite(energy)) {
            return 0;
        } else {
            return 1;
        }
    }

    public void setStdDev(double[][] o2, double[] coeff) {
        stdDev = new double[o2.length][o2[0].length];
        for (int i = 0; i < stdDev.length; i++) {
            for (int j = 0; j < stdDev[i].length; j++) {
                stdDev[i][j] = Math.sqrt(1.0 / (2.0 * o2[i][j] * coeff[i]));
            }
        }
    }
    
    
}
