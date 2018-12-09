package etomica.models.elasticconstantmappedavg;


import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataTensor;
import etomica.integrator.IntegratorMC;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeSumCrystal;
import etomica.lattice.crystal.*;
import etomica.math.function.Function;
import etomica.models.clathrates.molecularhma.Clathrateenergyandcv;
import etomica.molecule.IMoleculeList;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.normalmode.*;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.space3d.Space3D;
import etomica.space3d.Tensor3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Dimension;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;


public class variationwithstrain extends Simulation{

    public int[] nCells;
    public MeterPotentialEnergy meterPE;
    public Potential2SoftSpherical potentialLJ;
    public PrimitiveOrthorhombic primitiveOrthorhombic;
    public Basis basis;
    public BravaisLatticeCrystal lattice;
    public final Space space;
    public Box box;
//    public double[][] omega2;
 //   public double[][][] eigenvec;
    public final PotentialMasterCell potentialMaster;
     public IntegratorMC integrator;
    public SpeciesSpheresMono species;
    public double rCutLJ;
    public double harmonicE;


    public variationwithstrain(int i, double[] a0, double rCutLJ,int[] nCells, int[] nunitCells, Space space) {
         super(space);
         species = new SpeciesSpheresMono(this, space);
         species.setIsDynamic(true);
         addSpecies(species);
this.rCutLJ=rCutLJ;
        this.potentialLJ = potentialLJ;
        this.nCells = nCells.clone();
        this.space = space;
         Basis basis = new BasisCubicFcc();
         basis = new BasisBigCell(space, basis, nunitCells);
         primitiveOrthorhombic = new PrimitiveOrthorhombic(space, a0[0], a0[1], a0[2]);
//this.i=i;
         lattice = new BravaisLatticeCrystal(primitiveOrthorhombic, basis);
         Boundary boundary = new BoundaryDeformableLattice(primitiveOrthorhombic, nCells);
         box = this.makeBox(boundary);

         int numAtoms= 4*(nunitCells[0]*nunitCells[1]*nunitCells[2]);
   //      System.out.println(numAtoms);
         box.setNMolecules(species, numAtoms);
         potentialMaster = new PotentialMasterCell(this, rCutLJ, space);
         potentialMaster.setCellRange(2);

         if (rCutLJ > 0.49 * box.getBoundary().getBoxSize().getX(0)) {
             throw new RuntimeException("rc is too large");
         }

         potentialLJ = new P2LennardJones(space, 1, 1);
          AtomType leafType = species.getLeafType();
         potentialMaster.addPotential(potentialLJ, new AtomType[]{leafType, leafType});
         potentialMaster.getNbrCellManager(box).assignCellAll();
         System.out.println("Cell Density: " + numAtoms / boundary.volume());

         NormalModesPotential normalModesPotential=new NormalModesPotential(nCells, primitiveOrthorhombic, basis, potentialLJ, space);
         MCMoveHarmonicmodified mcMoveHarmonicmodified=new MCMoveHarmonicmodified();
         mcMoveHarmonicmodified.setWaveVectorCoefficients(normalModesPotential.getWaveVectorFactory().getCoefficients());
         mcMoveHarmonicmodified.setEigenVectors(normalModesPotential.getEigenvectors());
         mcMoveHarmonicmodified.setOmegaSquared(normalModesPotential.getOmegaSquared());
          CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitiveOrthorhombic, basis,space);
         coordinateDefinition.initializeCoordinates(nCells);
         CoordinateDefinition.BasisCell[] cells = coordinateDefinition.getBasisCells();
         mcMoveHarmonicmodified.setCoordinateDefinition(coordinateDefinition);
         mcMoveHarmonicmodified.setWaveVectors(normalModesPotential.getWaveVectorFactory().getWaveVectors());
         int coordinateDim=coordinateDefinition.getCoordinateDim();
         double[][] normalmodcoordreal = new double[normalModesPotential.getWaveVectorFactory().getWaveVectors().length][coordinateDim];
         normalmodcoordreal[0][0] = (i/100.0);
         int[] modnum= new int[12];
for(int hh=0;hh<modnum.length;hh++){modnum[hh]=hh;}

         System.out.println("normalmodcoordreal[0][0] "+normalmodcoordreal[0][0]);

                harmonicE= mcMoveHarmonicmodified.doTrial(normalmodcoordreal,cells,modnum,normalModesPotential.getOmegaSquared());

         System.out.println("harmonicE "+harmonicE);
     }



    public static void main(String[] args) {
        variationwithstrain.SimParams params = new variationwithstrain.SimParams();
        ParseArgs.doParseArgs(params, args);

         for(int j=0;j<20;j++){

         double[] a0={5+(j/20),5-(j/20),125/((5+(j/20))*(5-(j/20)))};
             System.out.println("this is for size: " + a0[0]+" "+ a0[1]+" "+ a0[2]);


            for(int i=0;i<100;i++){  //change normalmodecoord

            final variationwithstrain sim = new variationwithstrain(i, a0, params.rCutLJ, params.nCells,params.nunitCells, Space3D.getInstance());



            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
            meterPE.setBox(sim.box);
            System.out.println(meterPE.getData()+" ");}

        }
        //what is d above command doing
    //where is it taking normal mode coordinate
}
    public static class SimParams extends ParameterBase {
        public int[] nCells = {1,1,1};
        public int[] nunitCells = {4,4,4};
        public double rCutLJ = 1;

    }

}
