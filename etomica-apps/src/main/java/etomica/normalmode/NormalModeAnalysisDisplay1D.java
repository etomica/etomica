/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * 
 * 
 * @author Tai Boon Tan
 */
public class NormalModeAnalysisDisplay1D extends Simulation {

    private static final long serialVersionUID = 1L;
	private static final String APP_NAME = "1-D Harmonic Oscillator";

	public NormalModeAnalysisDisplay1D(Space space) {
        super(space);

        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        boundary = new BoundaryRectangularPeriodic(space, numAtoms / density);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        primitive = new PrimitiveCubic(space, 1.0 / density);
        nCells = new int[]{numAtoms};

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, space);
        coordinateDefinition.initializeCoordinates(nCells);

        nm = new NormalModes1DHR(boundary, numAtoms);
        nm.setTemperature(temperature);

        waveVectorFactory = nm.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);
        integrator = new IntegratorHarmonic(random, 0.01, temperature, box);
        integrator.setCoordinateDefinition(coordinateDefinition);
        integrator.setWaveVectors(waveVectorFactory.getWaveVectors());
        integrator.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        integrator.setOmegaSquared(nm.getOmegaSquared(), waveVectorFactory.getCoefficients());
        integrator.setEigenVectors(nm.getEigenvectors());
        integrator.setTemperature(temperature);


        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);

        getController().addAction(activityIntegrate);

    }
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
 
        //instantiate simulation
		Space sp = Space.getInstance(1);
        NormalModeAnalysisDisplay1D sim = new NormalModeAnalysisDisplay1D(sp);
        
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        simGraphic.makeAndDisplayFrame(APP_NAME);
        
	}
	
	public static class Applet extends javax.swing.JApplet{
		
        int numAtoms =30;
        double density = 1.0;
        double temperature = 0.1;
		
		public void init(){
	
		}
	}
	
	
	
	protected IntegratorHarmonic integrator;
	protected ActivityIntegrate activityIntegrate;
	protected Box box;
	protected Boundary boundary;
	protected Primitive primitive;
	protected Basis basis;
	protected int[] nCells;
	protected SpeciesSpheresMono species;
	protected NormalModes1DHR nm;
	protected CoordinateDefinitionLeaf coordinateDefinition;
	protected WaveVectorFactory waveVectorFactory;
	//Default:
	protected double density = 0.5;
	protected double temperature = 1.0;
	protected int numAtoms = 50;




	
}
