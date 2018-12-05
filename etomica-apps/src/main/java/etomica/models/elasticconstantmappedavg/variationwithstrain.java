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

     public variationwithstrain(double[] a0, double rCutLJ,int[] nCells, Space space) {
         super(space);
         species = new SpeciesSpheresMono(this, space);
         species.setIsDynamic(true);
         addSpecies(species);
this.rCutLJ=rCutLJ;
        this.potentialLJ = potentialLJ;
        this.nCells = nCells.clone();
        this.space = space;
         Basis basis = new BasisCubicFcc();
         primitiveOrthorhombic = new PrimitiveOrthorhombic(space, a0[0], a0[1], a0[2]);

         lattice = new BravaisLatticeCrystal(primitiveOrthorhombic, basis);
        int nSites = nCells[0] * nCells[1] * nCells[2];
        Boundary boundary = new BoundaryDeformableLattice(primitiveOrthorhombic, nCells);
        Box box = new Box(boundary, space);
        int numAtoms=nSites*4;
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
         System.out.println("Cell Density: " + nSites / boundary.volume());
    }



    public static void main(String[] args) {
        variationwithstrain.SimParams params = new variationwithstrain.SimParams();
        ParseArgs.doParseArgs(params, args);

         for(int j=0;j<20;j++){

         double[] a0={5+(j/20),5-(j/20),125/((5+(j/20))*(5-(j/20)))};

         final variationwithstrain sim = new variationwithstrain(a0, params.rCutLJ, params.nCells, Space3D.getInstance());


            for(int i=0;i<100;i++){  //change normalmodecoord


               NormalModesPotential normalModesPotential=new NormalModesPotential(params.nCells, sim.primitiveOrthorhombic, sim.basis, sim.potentialLJ, sim.space);

                MCMoveHarmonicmodified mcMoveHarmonicmodified=new MCMoveHarmonicmodified();
                int kDim = normalModesPotential.getWaveVectorFactory().getWaveVectors().length;
                double[] zero=new double [kDim];
                for(int z=0;z<kDim;z++) {zero[i] = 0;}
                mcMoveHarmonicmodified.setWaveVectorCoefficients(zero);
                mcMoveHarmonicmodified.setEigenVectors(normalModesPotential.getEigenvectors());
                mcMoveHarmonicmodified.setOmegaSquared(normalModesPotential.getOmegaSquared());

                int coordinateDim=3*4*params.nCells[0] * params.nCells[1] * params.nCells[2];
                double[][] normalmodcoordreal = new double[normalModesPotential.getWaveVectorFactory().getWaveVectors().length][coordinateDim];
                normalmodcoordreal[0][0] = (i/100.0);
                mcMoveHarmonicmodified.doTrial(normalmodcoordreal[0][0]);


            MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
            meterPE.setBox(sim.box);
            System.out.println(meterPE);}

        }
        //what is d above command doing
    //where is it taking normal mode coordinate
}
    public static class SimParams extends ParameterBase {
        //    	public String configFile = "config_from_paper_HHO_shiftedL_2_sI";
        public int[] nCells = {4,4,4};
        public double rCutLJ = 3;

    }

}
