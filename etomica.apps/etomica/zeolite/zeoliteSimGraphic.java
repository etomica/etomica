package etomica.zeolite;

import java.util.Iterator;
import java.util.LinkedList;

import javax.swing.JFrame;
import javax.swing.JPanel;

import etomica.graphics.DisplayPhase;
import etomica.action.Action;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.Atom;
import etomica.atom.AtomFilter;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.DeviceButton;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorPhase;
import etomica.integrator.IntervalActionAdapter;
import etomica.math.geometry.Plane;
import etomica.math.geometry.Polyhedron;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.simulation.SimulationContainer;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

public class zeoliteSimGraphic extends SimulationGraphic{
	public DisplayPhase display;
	public zeoliteSimGraphic(ZeoliteSimulation sim){
		super(sim);
		DeviceButtonSingle tSwitch = new DeviceButtonSingle(sim.getController());
		ZeoliteSimStart action = new ZeoliteSimStart(sim,this);
		tSwitch.setAction(action);
		tSwitch.setLabel("Start Simulation");
		super.add(tSwitch);
		
		display = super.getDisplayPhase(sim.phase);
		/*
		DeviceButton picture = new DeviceButton(sim.getController());
		TakePictures pic = new TakePictures(this);
		picture.setAction(pic);
		picture.setLabel("Picture Time");
		super.add(picture);
		*/
	}
	
}
