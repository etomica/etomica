package etomica.models.oneDHardRods;

import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.nbr.list.PotentialMasterList;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.normalmode.NormalModes;
import etomica.normalmode.NormalModes1DHR;
import etomica.normalmode.P2XOrder;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Null;


/**
 * Uses a Widom-like insertion of a mode to calculate a probability.
 * Uses a different box than the main simulation, to assume an extra mode & rod 
 * is added/
 * 
 * @author cribbin
 *
 */
public class MeterDifferentImageAdd extends DataSourceScalar {

    public int nInsert, counter;
    private MeterPotentialEnergy meterPE;
    private CoordinateDefinition cDef, simCDef;
    private int cDim, simCDim;
    private IVectorMutable[] waveVectors, simWaveVectors;
    private double[] simRealT, simImagT;
    private double[][] stdDev;
    protected double temperature;
    private double[] newU;
    private double[] wvCoeff, simWVCoeff;
    private double[][][] eigenVectors, simEigenVectors;

    protected final IRandom random;
    private IBox box;
    private int numAtoms;
    private Boundary bdry;
    private NormalModes nm;
    WaveVectorFactory waveVectorFactory;
    
    public MeterDifferentImageAdd(String string, /*IPotentialMaster potentialMaster,*/ 
            int numSimAtoms, double density, Simulation sim,
            Primitive simPrimitive, Basis simBasis, CoordinateDefinition simCD,
            NormalModes simNM, double temp){
        super(string, Null.DIMENSION);
        this.random = sim.getRandom();
        this.temperature = temp;
        
        simWaveVectors = simNM.getWaveVectorFactory().getWaveVectors();
        this.simCDef = simCD;
        simCDim = simCD.getCoordinateDim();
        simEigenVectors = simNM.getEigenvectors();
        simWVCoeff = simNM.getWaveVectorFactory().getCoefficients();
        simRealT = new double[simCDim];
        simImagT = new double[simCDim];
        
        numAtoms = numSimAtoms + 1;
        box = new Box(sim.getSpace());
        sim.addBox(box);
        box.setNMolecules(sim.getSpecies(0), numAtoms);
        bdry = new BoundaryRectangularPeriodic(sim.getSpace(), numAtoms/density);
        box.setBoundary(bdry);
        
        int[] nCells = new int[]{numAtoms};
        cDef = new CoordinateDefinitionLeaf(box, simPrimitive, 
                simBasis, sim.getSpace());
        cDef.initializeCoordinates(nCells);
        cDim = cDef.getCoordinateDim();

        nm = new NormalModes1DHR(box.getBoundary(), numAtoms);
        nm.setHarmonicFudge(1.0);
        nm.setTemperature(temperature);
        nm.getOmegaSquared();
        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        waveVectors = nm.getWaveVectorFactory().getWaveVectors();
        eigenVectors = nm.getEigenvectors();
        wvCoeff = nm.getWaveVectorFactory().getCoefficients();
        setStdDev(nm.getOmegaSquared(), wvCoeff);
        
        PotentialMasterList potentialMaster = new PotentialMasterList(sim, sim.getSpace());
        Potential2 potential = new P2HardSphere(sim.getSpace(), 1.0, true);
        potential = new P2XOrder(sim.getSpace(), (Potential2HardSpherical)potential);
        potential.setBox(box);
        potentialMaster.addPotential(potential, new IAtomType[] {((SpeciesSpheresMono)sim.getSpecies(0)).getLeafType(), ((SpeciesSpheresMono)sim.getSpecies(0)).getLeafType()});
        double neighborRange = 1.01/density;
        potentialMaster.setRange(neighborRange);
        //find neighbors now.  Don't hook up NeighborListManager since the
        //  neighbors won't change
        potentialMaster.getNeighborManager(box).reset();
        
        
        
        
        
        meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
    }
    
    public double getDataAsScalar() {
        
        IAtomList atomlist = box.getLeafList();
//        for (int i = 0; i < atomlist.getAtomCount(); i++){
//            System.out.println("start i " + atomlist.getAtom(i).getPosition().getX(0));
//        }

        BasisCell[] simCells = simCDef.getBasisCells();
        BasisCell[] cells = cDef.getBasisCells();
        BasisCell cell = simCells[0];
        //nan this makes it 1D
        newU = new double[cDim];
        for(int i=0; i < newU.length; i++){
            newU[i] = 0.0;
        }
        
        //Calculate normal mode coordinates of simulation system.
        double[] realCoord = new double[waveVectors.length];
        double[] imagCoord = new double[waveVectors.length];
        for (int wvcount = 0; wvcount < simWaveVectors.length; wvcount++){
            simCDef.calcT(simWaveVectors[wvcount], simRealT, simImagT);
            realCoord[wvcount] = 0.0;
            imagCoord[wvcount] = 0.0;
            for (int i = 0; i < simCDim; i++){
                for (int j = 0; j < simCDim; j++){
                    realCoord[wvcount] += simEigenVectors[wvcount][i][j] * simRealT[j];
                    imagCoord[wvcount] += simEigenVectors[wvcount][i][j] * simImagT[j];
                }
            }
            if(simWVCoeff[wvcount] != 1.0 ) {imagCoord[wvcount] = 0.0;}
        }
        
        //Create the last normal mode coordinate from the Gaussian distribution
        for (int j = 0; j < cDim; j++) {
            //We are adding 0.5, and this code lets us get it in the right slot.
            if(waveVectors.length == simWaveVectors.length){
                imagCoord[waveVectors.length - 1] = random.nextGaussian() * 
                    Math.sqrt(temperature) * stdDev[waveVectors.length - 1][j];
            } else {
                realCoord[waveVectors.length - 1] = random.nextGaussian() * 
                    Math.sqrt(temperature) * stdDev[waveVectors.length - 1][j];
            }
        }
            
        //Calculate the positions for the meter's system

        for (int iCell = 0; iCell < cells.length; iCell++){
            cell = cells[iCell];
            for (int j = 0; j < cDim; j++) {
                newU[j] = 0.0;
            }
            for (int wvcount = 0; wvcount < waveVectors.length; wvcount++){
                //Calculate the change in positions.
                double kR = waveVectors[wvcount].dot(cell.cellPosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for (int i = 0; i < cDim; i++){
                    for (int j = 0; j < cDim; j++){
                       newU[j] += wvCoeff[wvcount] * eigenVectors[wvcount][i][j] 
                            * 2.0 * (realCoord[wvcount] * coskR - imagCoord[wvcount] * sinkR);
                    }
                }
            }
            cDef.setToU(cells[iCell].molecules, newU);
        }
        
//        for (int i = 0; i < atomlist.getAtomCount(); i++){
//            System.out.println("end i " + atomlist.getAtom(i).getPosition().getX(0));
//        }
//        System.out.println(".");
        
        
        
        //Check for overlap & return based on it.
           if(Double.isInfinite(meterPE.getDataAsScalar())) {
//               System.out.println("0");
               return 0;
           } else {
//               System.out.println("1");
               return 1;
        }
    }

    private void setStdDev(double[][] o2, double[] coeff) {
        stdDev = new double[o2.length][o2[0].length];
        for (int i = 0; i < stdDev.length; i++) {
            for (int j = 0; j < stdDev[i].length; j++) {
                stdDev[i][j] = Math.sqrt(1.0 / (2.0 * o2[i][j] * coeff[i]));
            }
        }
    }
    
    
}
