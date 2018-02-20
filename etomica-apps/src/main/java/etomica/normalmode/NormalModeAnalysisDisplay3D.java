/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;

/**
 * 
 * 
 * @author Tai Boon Tan
 */
public class NormalModeAnalysisDisplay3D extends Simulation {

    private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "3-D Harmonic Oscillator";
    protected IntegratorHarmonic integrator;
    protected ActivityIntegrate activityIntegrate;
    protected Box box;
    protected Boundary boundary;
    protected Primitive primitive;
    protected Basis basis;
    protected int[] nCells;
    protected SpeciesSpheresMono species;
    protected NormalModes3D nm;
    protected WaveVectorFactory waveVectorFactory;
    protected CoordinateDefinitionLeaf coordinateDefinition;
    protected MeterPotentialEnergy meterPE;
    protected Space space;
    protected double L;
    protected int n = 3;
    protected P2SoftSphericalTruncatedShifted pTruncated;
    protected double truncationRadius;
    protected double density = 1.256;
    protected double temperature = 0.0001;
    protected double latticeEnergy;

    public NormalModeAnalysisDisplay3D(Space _space) {
        super(_space);
        this.space = _space;

        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        box = this.makeBox();
        box.setNMolecules(species, 4 * n * n * n);

        L = Math.pow(4.0 / density, 1.0 / 3.0);
        primitive = new PrimitiveCubic(space, L);

        nCells = new int[]{n, n, n};
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        basis = new BasisCubicFcc();
        box.setBoundary(boundary);

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, 12);
        truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
        pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        AtomType sphereType = species.getLeafType();

        PotentialMasterMonatomic potentialMaster = new PotentialMasterMonatomic(this);
        potentialMaster.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});
        meterPE = new MeterPotentialEnergy(potentialMaster, box);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);

        latticeEnergy = meterPE.getDataAsScalar();
        //String fileName = "CB_FCC_n12_T01_Mode01";
        nm = new NormalModes3D(space, primitive, basis);
        nm.setTemperature(temperature);
        nm.setNCellNum(n);


        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        integrator = new IntegratorHarmonic(random, 0.0001, temperature, box);

        integrator.setOmegaSquared(nm.getOmegaSquared(), waveVectorFactory.getCoefficients());
        integrator.setEigenVectors(nm.getEigenvectors());
        integrator.setWaveVectors(waveVectorFactory.getWaveVectors());
        integrator.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        integrator.setCoordinateDefinition(coordinateDefinition);
        integrator.setTemperature(temperature);
        //integrator.setOneWV(true);
        //integrator.setWaveVectorNum(0);
        //integrator.setOneEVal(true);
        //integrator.setEValNum(4);


        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(0);

        getController().addAction(activityIntegrate);

    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        //instantiate simulation
        NormalModeAnalysisDisplay3D sim = new NormalModeAnalysisDisplay3D(Space.getInstance(3));

        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
        simGraphic.makeAndDisplayFrame(APP_NAME);

    }

	public int[] getNCells() {
		return nCells;
	}

	public void setNCells(int[] cells) {
		nCells = cells;
	}

    public int getN() {
		return n;
	}

	public void setN(int n) {
		this.n = n;
	}

	public static class Applet extends javax.swing.JApplet{

		public void init(){

		}
	}
	
}
