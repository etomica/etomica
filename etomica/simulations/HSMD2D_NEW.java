package etomica.simulations;

import java.awt.Color;
import java.awt.GridLayout;

import etomica.Action;
import etomica.Atom;
import etomica.Controller;
import etomica.Default;
import etomica.Integrator;
import etomica.IntegratorIntervalEvent;
import etomica.IntegratorIntervalListener;
import etomica.Modifier;
import etomica.Phase;
import etomica.Simulation;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseQuench;
import etomica.data.DataSourceScalar;
import etomica.graphics.ColorScheme;
import etomica.graphics.DeviceButton;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceToggleButton;
import etomica.graphics.DeviceTrioControllerButton;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.log.LoggerAbstract;
import etomica.modifier.ModifierBoolean;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential2;
import etomica.space2d.Space2D;
import etomica.units.Dimension;


/**
 * @author kofke
 *
 * Non-equilibrium work calculation of the pulling of an atom through a fluid.
 * Exponential work for process should be zero.
 */
public class HSMD2D_NEW extends SimulationGraphic  {

	/**
	 * @author kofke
	 *
	 * To change this generated comment edit the template variable "typecomment":
	 * Window>Preferences>Java>Templates.
	 * To enable and disable the creation of type comments go to
	 * Window>Preferences>Java>Code Generation.
	 */
	public class MyBoolModifier implements ModifierBoolean {

		/**
		 * @see etomica.modifier.ModifierBoolean#setBoolean(boolean)
		 */
		public void setBoolean(boolean b) {auto = b;}

		/**
		 * @see etomica.modifier.ModifierBoolean#getBoolean()
		 */
		public boolean getBoolean() {
			return auto;
		}

	}
	public class MyAction extends Action implements IntegratorIntervalListener {
		boolean moving = false;
		boolean movingDown = true;
		double ke0, ke1;
		int waitCount;
		PhaseQuench quencher = new PhaseQuench(phase);
        public int getPriority() {return 300;}
		public void actionPerformed() {
			if(!moving) {//start motion
				first.coord.momentum().E((movingDown?+1:-1)*speed*mass);
				ke0 = first.coord.momentum().squared()/2.0/mass;
				integrator.reset();
				moving = true;
			} else {//stop motion
				doStop();
			}
		}
		private void doStop() {
			ke1 = first.coord.momentum().squared()/2.0/mass;
			System.out.println("ke0, ke1, w: "+ke0+" "+ke1+" "+(ke1-ke0));
			wMeter.updateSums();
			eMeter.updateSums();
			eHist.doUpdate();
			wHist.doUpdate();
			eHist.repaint();
			wHist.repaint();
			first.coord.momentum().E(0.0);
			phase.getAgent(species).coord.accelerateTo(space.makeVector());
			quencher.actionPerformed();
			moving = false;
			movingDown = !movingDown;
			if(auto) waitCount = 10;			
		}
		public void intervalAction(IntegratorIntervalEvent evt) {
			if(auto) {
				if(waitCount > 0) waitCount--;
				else if(!moving) {actionPerformed(); return;} 
			}
			if(!moving) return;
			double x = first.coord.position().x(0);
			double y = first.coord.position().x(1);
			double d = phase.boundary().dimensions().x(0);
			if((!movingDown && (x < 0.05*d || y < 0.05*d)) || (movingDown && (x > 0.95*d || y > 0.95*d))) doStop();
		}
	}
	public class MyColorScheme extends ColorScheme {
		public Color atomColor(Atom a) {
			return (a == first) ? java.awt.Color.red : java.awt.Color.black;
		}
	}
	public IntegratorHard integrator;
	public SpeciesSpheresMono species;
	public Phase phase;
	public Potential2 potential;
	public Controller controller;
	public DisplayPhase display;
	public DisplayPlot eHist, wHist;
	public boolean auto;
	Atom first;
	double speed;
	double mass = 1e10;
	MyMeter wMeter, eMeter;

	public HSMD2D_NEW() {
		super(Space2D.getInstance());
		etomica.units.BaseUnit.Length.Sim.TO_PIXELS = 8;
		Simulation.instance = this;
		integrator = new IntegratorHard(this);
		integrator.setSleepPeriod(2);
		integrator.setInterval(1);
//		integrator.setIsothermal(true);
		species = new SpeciesSpheresMono(this);
		species.setNMolecules(60);
		phase = new Phase(this);
		potential = new P2HardSphere();
		potential.setSpecies(species,species);
		controller = new Controller(this);
		DeviceTrioControllerButton controllerButtons = new DeviceTrioControllerButton(this, controller);
		display = new DisplayPhase(this);
		display.setColorScheme(new MyColorScheme());
		MyAction action = new MyAction();
		DeviceButton button = new DeviceButton(action);
		integrator.addListener(action);
		DeviceSlider speedSlider = new DeviceSlider(new MyModifier());
		speedSlider.setMinimum(0.0);
		speedSlider.setMaximum(200.0);
		speedSlider.setValue(10.0);
		wMeter = new MyMeter(this,false,action);
		eMeter = new MyMeter(this,true,action);
		eMeter.setHistogramming(true);
		wMeter.setHistogramming(true);
		eHist = new DisplayPlot();
		wHist = new DisplayPlot();
//		panel().setBackground(java.awt.Color.yellow);
//		elementCoordinator.go();
		first = species.getAgent(phase).firstMolecule();
		first.coord.setMass(mass);
		etomica.utility.HistogramSimple eHistogram = eMeter.getHistogram();
		etomica.utility.HistogramSimple wHistogram = wMeter.getHistogram();
		eHist.setDataSources(eHistogram);
		wHist.setDataSources(wHistogram);
		eHistogram.setXLabel("exp(-W/kT)");
		wHistogram.setXLabel("W/kT");
		eHist.setXLabel("exp(-W/kT)");
		wHist.setXLabel("W/kT");
		
		display.setSize(200,200);
		eHist.setSize(300,200);
		wHist.setSize(300,200);
		eHist.setUpdateInterval(1000);
		wHist.setUpdateInterval(1000);
		eHist.getPlot().setXLabel("exp(-W/kT)");
		wHist.getPlot().setXLabel("W/kT");
		DeviceToggleButton auto = new DeviceToggleButton(this,new MyBoolModifier(),"Auto","Manual");
		
		panel().remove(panel().devicePanel);
		javax.swing.JPanel myDevicePanel = new javax.swing.JPanel();
		myDevicePanel.setLayout(new GridLayout(1,3));
		javax.swing.JPanel deviceSubPanel = new javax.swing.JPanel();
		deviceSubPanel.setLayout(new GridLayout(2,1));
		deviceSubPanel.add(button.graphic());
		deviceSubPanel.add(speedSlider.graphic());
		myDevicePanel.add(controllerButtons.graphic());
		myDevicePanel.add(deviceSubPanel);
		myDevicePanel.add(auto.graphic());
		panel().add(myDevicePanel);
		
		panel().remove(panel().displayPanel);
		javax.swing.JPanel myPanel = new javax.swing.JPanel();
		myPanel.setLayout(new GridLayout(1,2));
		panel().add(myPanel);
		myPanel.add(display.graphic());
		javax.swing.JPanel histPanel = new javax.swing.JPanel();
		histPanel.setLayout(new GridLayout(2,1));
		histPanel.add(eHist.graphic());
		histPanel.add(wHist.graphic());
		myPanel.add(histPanel);
		
		Logger logger = new Logger(this, integrator);
		logger.setUpdateInterval(Default.BLOCK_SIZE);
		logger.setAppending(false);
		logger.setCloseFileEachTime(true);

	}
	
	class MyModifier implements Modifier {
		public void setValue(double v) {speed = (v+1)/10.;}
		public double getValue() {return speed*10-1;}
		public etomica.units.Dimension getDimension() {return etomica.units.Dimension.UNDEFINED;}
		public String getLabel() {return "speed";}
    }
	
	class MyMeter extends DataSourceScalar {

		boolean doExp;
		MyAction action;
		/**
		 * Constructor for MyMeter.
		 * @param parent
		 */
		public MyMeter(SimulationElement parent, boolean doExp, MyAction action) {
			super(parent);
			this.doExp = doExp;
			this.action = action;
		}

		/**
		 * @see etomica.data.DataSourceScalar#getData()
		 */
		public double getDataAsScalar(Phase p) {
			double w = -(action.ke1 - action.ke0)/HSMD2D_NEW.this.integrator.getTemperature();
			if(!doExp) System.out.println(w);
			return doExp ? ((w>1000)?1e-10:Math.exp(-w)) : w;
		}

		/**
		 * @see etomica.data.Meter#getDimension()
		 */
		public Dimension getDimension() {
			return Dimension.NULL;
		}
	}//end MyMeter

	public class Logger extends LoggerAbstract {

		/**
		 * Constructor for Logger.
		 * @param sim
		 * @param integrator
		 */
		public Logger(Simulation sim, Integrator integrator) {
			super(sim, integrator);
		}

		protected void write() throws java.io.IOException {
			int n = wMeter.getHistogram().getHistoryLength();
			double[] wVal = wMeter.getHistogram().getData(null);
			double[] eVal = eMeter.getHistogram().getData(null);
			double[] xw = wMeter.getHistogram().getX();
			double[] xe = eMeter.getHistogram().getX();
			for(int i=0; i<n; i++) {
				fileWriter.write(Double.toString(xw[i])+"   "+
								 Double.toString(wVal[i])+"   "+
								 Double.toString(xe[i])+"   "+
								 Double.toString(eVal[i])+"\n");
			}
			System.out.println("writing");
		}

	}
	
	
    
	/**
	 * Demonstrates how this class is implemented.
	 */
	public static void main(String[] args) {
		HSMD2D_NEW sim = new HSMD2D_NEW();
		SimulationGraphic.makeAndDisplayFrame(sim);
	//	sim.controller.start();
	}//end of main



}
