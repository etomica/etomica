//Source file generated by Etomica

package etomica.simulation.prototypes;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPhase;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicFcc;
import etomica.modifier.Modifier;
import etomica.phase.Phase;
import etomica.potential.P2SquareWell;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Dimension;
import etomica.units.Length;

//remember to set up Space3D.CoordinatePair.reset if experiencing 
//problems with this simulation hanging

public class SWMD3D extends Simulation {

	public class MyModifier implements Modifier {

        public Dimension getDimension() {
            return Length.DIMENSION;
        }

		/**
		 * @see etomica.modifier.Modifier#setValue(double)
		 */
		public void setValue(double d) {
			potential.setCoreDiameter(d);
			((AtomTypeSphere)species.getMoleculeType()).setDiameter(d);
		}

		/**
		 * @see etomica.modifier.Modifier#getValue()
		 */
		public double getValue() {
			return potential.getCoreDiameter();
		}
        
        public String getLabel() {
            return "diameter";
        }

	}
	
  public SWMD3D() {
	super(Space3D.getInstance());
//	defaults.makeLJDefaults();
    
    integrator = new IntegratorHard(this);
//  integrator.addIntervalListener(((PotentialMasterNbr)potentialMaster).getNeighborManager());
    integrator.setTimeStep(0.01);
    integrator.setIsothermal(true);
    integrator.setTemperature(300);
    ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
    getController().addAction(activityIntegrate);


    phase = new Phase(this);
    potential  = new etomica.potential.P2SquareWell(this);
    potential.setLambda(1.6);

    species  = new etomica.species.SpeciesSpheresMono(this);
    getSpeciesRoot().addSpecies(species);
    phase.getAgent(species).setNMolecules(108);

	
//	DeviceSlider tControl = new DeviceSlider(integrator, "temperature");
//	DeviceSlider sigmaControl = new DeviceSlider(new MyModifier());
//	DeviceSlider lambdaControl = new DeviceSlider(potential0, "lambda");
//	tControl.setLabel("Temperature (K)");
//	sigmaControl.setLabel("Atom size (Angstroms)");
//	tControl.setShowValues(true);
//	tControl.setShowBorder(true);
//	tControl.setMinimum(100);
//	tControl.setMaximum(700);
//	sigmaControl.setShowValues(true);
//	sigmaControl.setShowBorder(true);
//	sigmaControl.setPrecision(2);
//	sigmaControl.setMinimum(0.0);
//	sigmaControl.setMaximum(3.0);
//	lambdaControl.setShowValues(true);
//	lambdaControl.setShowBorder(true);
//	lambdaControl.setPrecision(2);
//	lambdaControl.setMinimum(1.1);
//	lambdaControl.setMaximum(2.1);
//	lambdaControl.setValue(1.4);
//	lambdaControl.setNMajor(5);


//	mediator().go();
    this.potentialMaster.addPotential(potential,new Species[]{species,species});

    integrator.setPhase(phase);
    integrator.addListener(new IntervalActionAdapter(new PhaseImposePbc(phase)));

//	DeviceNSelector nControl = new DeviceNSelector(speciesSpheres0.getAgent(phase0));
//	nControl.setMaximum(108);
	phase.setDensity(0.0405);
    ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc());
    configuration.initializeCoordinates(phase);
  }

  private static final long serialVersionUID = 1L;
  public IntegratorHard integrator;
  public SpeciesSpheresMono species;
  public Phase phase;
  public P2SquareWell potential;
  public Controller controller;
  public DisplayPhase display;


  
  public static class MyColorScheme extends ColorScheme {
      public MyColorScheme(AtomLeaf redAtom) {
          atom = redAtom;
      }
	  public java.awt.Color getAtomColor(AtomLeaf a) {
		  return (a == atom) ? java.awt.Color.red : java.awt.Color.yellow;
	  }
      private static final long serialVersionUID = 1L;
      private AtomLeaf atom;
  }

}//end of class
