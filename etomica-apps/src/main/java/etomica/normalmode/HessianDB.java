/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Tensor3D;
import etomica.species.SpeciesSpheresMono;

import java.io.FileWriter;
import java.io.IOException;

/**
 * MC simulation of FCC soft-sphere model in 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 * 
 * @author Tai Boon Tan
 */
public class HessianDB extends Simulation {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
    final Tensor3D identity = new Tensor3D(new double[][]{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
    final Tensor3D tensor = new Tensor3D();
    public Box box;
    public Boundary boundary;
    public Primitive primitive, primitiveUnitCell;
    public Basis basis, basisFCC;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public int[] nCells;
    public MCMoveHarmonic move;
    public CoordinateDefinition coordinateDefinition;
    public PotentialMasterMonatomic potentialMaster;
    public NormalModes nm;
    protected double latticeEnergy;

    public HessianDB(Space _space, int numAtoms, double density, double temperature, int exponent, String filename) {
        super(_space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        /*
         * Creating new basis
         */

        potentialMaster = new PotentialMasterMonatomic(this);
        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        integrator = new IntegratorMC(this, potentialMaster, box);
        box.setNMolecules(species, numAtoms);

        primitive = new PrimitiveCubic(space, n * L);
        nCells = new int[]{n, n, n};
        basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, 12);
        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
        System.out.println("radius: " + truncationRadius);
        P2SoftSphericalTruncated pTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        Vector pos1 = space.makeVector();
        Vector pos2 = space.makeVector();
        Vector r = space.makeVector();


        try {
            FileWriter fileWriterH = new FileWriter(filename + ".h");
            FileWriter fileWriterVal = new FileWriter(filename + ".val");
            FileWriter fileWriterVec = new FileWriter(filename + ".vec");
            FileWriter fileWriterK = new FileWriter(filename + ".k");

            WaveVectorFactory wv = new WaveVectorFactorySimple(primitive, space);
            wv.makeWaveVectors(box);

            int rdim = numAtoms * space.D();

            double[][] array = new double[rdim][rdim];

            for (int atomN1 = 0; atomN1 < numAtoms; atomN1++) {
                for (int atomN2 = 0; atomN2 < numAtoms; atomN2++) {
                    if (atomN2 == atomN1) continue;

                    IAtom atom1 = box.getLeafList().get(atomN1);
                    IAtom atom2 = box.getLeafList().get(atomN2);

                    pos1 = atom1.getPosition();
                    pos2 = atom2.getPosition();

                    r.Ev1Mv2(pos2, pos1);

                    box.getBoundary().nearestImage(r); // get the nearest image

                    double[][] der2 = new double[space.D()][space.D()];
                    derivative2nd(r, pTruncated).assignTo(der2);
                    for (int i = 0; i < space.D(); i++) {
                        for (int j = 0; j < space.D(); j++) {

                            array[atomN1 * space.D() + i][atomN2 * space.D() + j] = der2[i][j];
                        }
                    }

                }

            }

            // self-term
            for (int atomN1 = 0; atomN1 < numAtoms; atomN1++) {
                for (int atomN2 = 0; atomN2 < numAtoms; atomN2++) {
                    if (atomN1 == atomN2)
                        continue; // we might double sum the elements in array[a][a] if we don't skip the pair
                    for (int i = 0; i < space.D(); i++) {
                        for (int j = 0; j < space.D(); j++) {

                            array[atomN1 * space.D() + i][atomN1 * space.D() + j] -= array[atomN1 * space.D() + i][atomN2 * space.D() + j];
                        }
                    }
                }
            }

            /*
             * impose symmetry on the matrix
             * Jama would generate a non-orthogonal eigenvectors if the last-digit
             *  in value in the matrix is different (numerical precision problem).
             */
            double[][] arrayjjp = array.clone();
            double[][] arrayjpj = array.clone();

            for (int i = 0; i < array.length; i++) {
                for (int j = i + 1; j < array[0].length; j++) {
                    double ave = 0.5 * (arrayjjp[i][j] + arrayjpj[j][i]);
                    array[i][j] = ave;
                    array[j][i] = ave;
                }
            }

            for (int i = 0; i < space.D() * numAtoms; i++) {
                for (int j = 0; j < space.D() * numAtoms; j++) {
                    fileWriterH.write(array[i][j] + " ");
                }
                fileWriterH.write("\n");
            }


            Matrix matrix = new Matrix(array);
            EigenvalueDecomposition ed = new EigenvalueDecomposition(matrix);

            double[] eVals = ed.getRealEigenvalues();
            double[][] eVecs = ed.getV().getArray();
            double[] kCoefficients = wv.getCoefficients();


            // output .k file
            for (int i = 0; i < kCoefficients.length; i++) {
                fileWriterK.write(Double.toString(kCoefficients[i]));
                for (int j = 0; j < wv.getWaveVectors()[i].getD(); j++) {
                    fileWriterK.write(" " + wv.getWaveVectors()[i].getX(j));

                }
                fileWriterK.write("\n");
            }
            // output .val file
            for (int ival = 0; ival < eVals.length; ival++) {
                if (eVals[ival] < 1E-12) {
                    fileWriterVal.write("0.0 ");
                } else {
                    fileWriterVal.write(1 / eVals[ival] + " ");
                }
            }


            // output .vec file
            for (int ivec = 0; ivec < rdim; ivec++) {
                for (int jvec = 0; jvec < rdim; jvec++) {
                    if (Math.abs(eVecs[jvec][ivec]) < 1e-15) {
                        fileWriterVec.write("0.0 ");
                    } else {
                        fileWriterVec.write(eVecs[jvec][ivec] + " ");
                    }
                }
                fileWriterVec.write("\n");
            }

            fileWriterH.close();
            fileWriterVal.close();
            fileWriterVec.close();
            fileWriterK.close();


        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA =500;
        double density = 1.256;
        double temperature = 1.0;
        int exponent = 12 ;
        String filename = "inputSSDB"+nA;//+"_d0962";
        // construct simulation
        HessianDB sim = new HessianDB(Space.getInstance(D), nA, density, temperature, exponent, filename);


        // SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
       // simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
       // simGraphic.makeAndDisplayFrame("Test");


//        MeterWorkHarmonicPhaseSpace meterWorkHarmonic = new MeterWorkHarmonicPhaseSpace(sim.move, sim.potentialMaster);
//        meterWorkHarmonic.setTemperature(temperature);
//        meterWorkHarmonic.setLatticeEnergy(sim.latticeEnergy);
//
//        AccumulatorAverageFixed dataAverage = new AccumulatorAverageFixed();
//        DataPump pump = new DataPump(meterWorkHarmonic, dataAverage);
//        IntegratorListenerAction pumpActionListener = new IntegratorListenerAction(pump);
//        pumpActionListener.setInterval(1);
//        sim.integrator.getEventManager().addListener(pumpActionListener);
//
//        sim.activityIntegrate.setMaxSteps(1000);
//        sim.getController().actionPerformed();
    }

    public Tensor3D derivative2nd(Vector r, Potential2SoftSpherical potential) {

        tensor.Ev1v2(r, r);
        double r2 = r.squared();
        double dW = potential.du(r2);
        double d2W = potential.d2u(r2);
        tensor.TE(1.0 / (r2 * r2) * (dW - d2W));
        tensor.PEa1Tt1(-dW / r2, identity);

        return tensor;
    }
    
}
