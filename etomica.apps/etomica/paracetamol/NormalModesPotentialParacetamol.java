package etomica.paracetamol;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.box.Box;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeSumCrystal.DataGroupLSC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.NormalModes;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.IVector;
import etomica.units.Dimension;
import etomica.util.Arrays;
import etomica.util.Function;
import etomica.util.FunctionGeneral;
import etomica.util.RandomNumberGenerator;
import etomica.util.numerical.FiniteDifferenceDerivative;

/**
 * Uses analysis of 2nd derivatives to compute the normal modes for a Bravais lattice with a basis, 
 * occupied by atoms that interact with a simple, spherically-symmetric soft potential.
 */

public class NormalModesPotentialParacetamol implements NormalModes {

    public NormalModesPotentialParacetamol(int[] nCells, Primitive primitive, Basis basis) {
        
        harmonicFudge = 1.0;
        needToCalculateModes = true;
        
        lattice = new BravaisLatticeCrystal(primitive, basis);
        
        int nSites = nCells[0]*nCells[1]*nCells[2];
        Boundary boundary = new BoundaryDeformableLattice(primitive, new RandomNumberGenerator(), nCells);
        
        box = new Box(boundary);

        System.out.println("Cell Density: "+nSites/boundary.volume());
        System.out.println("Site Density: "+nSites/boundary.volume()*basis.getScaledCoordinates().length);
        kFactory = new WaveVectorFactorySimple(primitive);
        kFactory.makeWaveVectors(box);
        System.out.println("Number of wave vectors: "+kFactory.getWaveVectors().length);
        
        double sum = 0.0;
        for(int i=0; i<kFactory.getWaveVectors().length; i++) {
            sum += kFactory.getCoefficients()[i];
        }
        System.out.println("nCells: "+Arrays.toString(nCells));
        System.out.println("Number of wave vectors represented: "+2.0*sum);
    }
    
    public void calculateModes() {
    	
        int kDim = kFactory.getWaveVectors().length;
        int spaceDim = lattice.getSpace().D();
        int basisDim = lattice.getBasis().getScaledCoordinates().length;
        int eDim = basisDim * spaceDim;
        omega2 = new double[kDim][eDim];
        eigenvectors = new double[kDim][eDim][eDim];
        
        
        //this function returns phi_{alpha,beta}, as defined in Dove Eq. 6.15
        FunctionGeneral function = new FunctionGeneral() {
            public Data f(Object obj) {
            	//Object is passed in as a set of generalized coordinates within a cell
                
            	DataDoubleArray u = (DataDoubleArray)obj;
            	
            	int[] d2 = new int[coordinateDefinitionParacetamol.getCoordinateDim()];
            	
            	////////////////////////////////////////////////////////////
//          	Vector3D r = (Vector3D)obj;
//              tensor.x.Ev1v2(r, r);
//              double r2 = r.squared();
//              double dW = potential.du(r2);
//              double d2W = potential.d2u(r2);
//              tensor.TE(1.0/(r2*r2)*(dW - d2W));
//              tensor.x.PEa1Tt1(-dW/r2,identity);
//      		  return tensor;
             
              //////////////////////////////////////////////////////////////
            	NumericalDerivativeEnergyParacetamol dW = new NumericalDerivativeEnergyParacetamol(box, potentialMaster);
            	dW.setCoordinateDefinition(coordinateDefinitionParacetamol);
            	
            	for (int i=0; i<coordinateDefinitionParacetamol.getCoordinateDim(); i++){
            		
            		/*
            		 * assignment of the first derivative
            		 */
            		d2[i]++;
            		double[] columnArrayi = new double[coordinateDefinitionParacetamol.getCoordinateDim()];
            		
            		for (int j=0; j<coordinateDefinitionParacetamol.getCoordinateDim(); j++){
            			/*
                		 * assignment of the second derivative
                		 */
            			d2[j]++;
            			
		            	FiniteDifferenceDerivative d2W = new FiniteDifferenceDerivative(dW);
		        		d2W.setH(0.00001);
		        		d2W.setHOptimizer(false);
		        		d2W.setNtab(10);
		        		
		        		columnArrayi[j] = d2W.df(d2, u.getData());
		        		d2[j]--;
		        	}
            		tensorDoubleArray.assignColumnFrom(i, columnArrayi);
            		d2[i]--;
            	}
            	return tensorDoubleArray;
                
            }
            
            public IDataInfo getDataInfo() {
                return dataInfo;
            }
            
            final DataDoubleArray tensorDoubleArray = new DataDoubleArray(new int[]{coordinateDefinitionParacetamol.getCoordinateDim(),
					   coordinateDefinitionParacetamol.getCoordinateDim()});;
            //final DataTensor tensor = new DataTensor(Space3D.getInstance());
            final DataInfo dataInfo = new DataDoubleArray.DataInfoDoubleArray("",Dimension.MIXED.MIXED, new int[]{coordinateDefinitionParacetamol.getCoordinateDim(),
					   coordinateDefinitionParacetamol.getCoordinateDim()});
            //final Tensor3D identity = new Tensor3D(new double[] {1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0});
       
        
        };
        
        LatticeSumCrystalParacetamol summer = new LatticeSumCrystalParacetamol(lattice);
        summer.setMaxLatticeShell(maxLatticeShell);
        IVector kVector = lattice.getSpace().makeVector();
        
        //calculation of self term
        kVector.E(0.0);
        summer.setK(kVector);
        System.out.println("\n k:"+kVector.toString()+"   in NormalModesPotentialParacetamol");
        summer.setCoordinateDefinitionParacetamol(coordinateDefinitionParacetamol);
        DataGroupLSC sum0 = (DataGroupLSC)summer.calculateSum(function);
        Function chopper = new Function.Chop(1e-9);
        sum0.map(chopper);
//        System.out.println(sum0.toString());
        double[][] array = new double[eDim][eDim];
        for(int k=0; k<kDim; k++) {
            kVector.E(kFactory.getWaveVectors()[k]);
            summer.setK(kVector);
            System.out.println("k:"+kVector.toString());
            DataGroupLSC sum = (DataGroupLSC)summer.calculateSum(function);
            sum.map(chopper);
            for(int j=0; j<basisDim; j++) {
                for(int jp=0; jp<basisDim; jp++) {
                    sum.getDataReal(j,j).ME(sum0.getDataReal(j, jp));
                }
            }
            for(int j=0; j<basisDim; j++) {
                for(int jp=0; jp<basisDim; jp++) {
                	
                	
                	/////////////////////////
                	DataDoubleArray tensor = ((DataDoubleArray)sum.getDataReal(j, jp));
                    //Tensor tensor = ((DataTensor)sum.getDataReal(j,jp)).x;
                    
                    
                    for(int alpha=0; alpha<spaceDim; alpha++) {
                        for(int beta=0; beta<spaceDim; beta++) {
                        	
                        	//////////////////////////
                            array[spaceDim*j+alpha][spaceDim*jp+beta] = tensor.getValue(new int[]{alpha, beta});
                        }
                    }
                }
            }
            Matrix matrix = new Matrix(array);
            EigenvalueDecomposition ed = matrix.eig();
            double[] eVals = ed.getRealEigenvalues();
            double[][] eVecs = ed.getV().getArray();
            System.out.println("Real eigenvalues: " + Arrays.toString(eVals));
//            System.out.println("Imag eigenvalues: " + Arrays.toString(ed.getImagEigenvalues()));
            
            for(int j=0; j<eDim; j++) {
                omega2[k][j] = eVals[j];
                for(int m=0; m<eDim; m++) {
                    eigenvectors[k][j][m] = eVecs[m][j];//need to check if indexes are right
                }
            }
//            System.out.println(sum[0].toString());
//            System.out.println();
//            System.out.println(sum[1].toString());
        }
        needToCalculateModes = false;
    }
 

    public double[][][] getEigenvectors(Box box) {
        if(needToCalculateModes) {
            calculateModes();
        }
        return eigenvectors;
    }

    public double[][] getOmegaSquared(Box box) {
        if(needToCalculateModes) {
            calculateModes();
        }
        return omega2;
    }

    public WaveVectorFactory getWaveVectorFactory() {
        return kFactory;
    }

    public int getMaxLatticeShell() {
        return maxLatticeShell;
    }

    public void setMaxLatticeShell(int maxLatticeShell) {
        this.maxLatticeShell = maxLatticeShell;
    }

    public void setHarmonicFudge(double newHarmonicFudge) {
        needToCalculateModes = true;
        harmonicFudge = newHarmonicFudge;
    }
    
    public void setTemperature(double newTemperature) {
        needToCalculateModes = true;
        temperature = newTemperature;
    }

    public void setCoordinateDefinition(CoordinateDefinitionParacetamol coordinateDefinitionParacetamol ){
    	this.coordinateDefinitionParacetamol = coordinateDefinitionParacetamol;
    }
    
    public CoordinateDefinitionParacetamol getCoordinateDefinition(){
    	return coordinateDefinitionParacetamol;
    }
    
	public void setPotentialMaster(PotentialMaster potentialMaster) {
		this.potentialMaster = potentialMaster;
	}
    
	public Box getBox() {
		return box;
	}

	public void setBox(Box box) {
		this.box = box;
	}
	
    private CoordinateDefinitionParacetamol coordinateDefinitionParacetamol;
    private PotentialMaster potentialMaster;
    private Box box;
    private final BravaisLatticeCrystal lattice;
    private Potential2SoftSpherical potential;
    private WaveVectorFactory kFactory;
    private int maxLatticeShell;
    private double harmonicFudge;
    private double temperature;
    private double[][] omega2;
    private double[][][] eigenvectors;
    private boolean needToCalculateModes;
    
	public PotentialMaster getPotentialMaster() {
		return potentialMaster;
	}





}
