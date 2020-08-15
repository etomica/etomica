/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.zeolite;

import etomica.action.IAction;
import etomica.action.activity.ControllerEvent;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMD;
import etomica.space.Vector;
import etomica.util.IEvent;
import etomica.util.IListener;

import java.io.FileWriter;
import java.io.IOException;

/* =====SUMMARY======
 * At each 'writeInterval', which corresponds to a certain number of simulation steps,
 * every atom's absolute distance traveled is written to a file of the name 'fileName'.
 * Subsequent 'writeInverval' outputs append the file to create repeating blocks of the 
 * simulation's atoms with columns of X, Y, and Z.
 * ==================
 * 
 * The class is ordered as follows.  After the atoms have moved, but before periodic
 * boundary conditions are imposed, the main class's 'intervalAction' method gets atom
 * coordinates and stores them in 'atomOldCoord'.  The 'intervalCount', initially set to
 * 'writeInterval', is predecremented each simulation step.  Each time it does NOT equal
 * zero, the class continues and 'intervalAction' of the subclass 'AfterPBC' is called.  
 * The method fills 'atomPBIarray' using the pre-periodic boundary* 'atomOldCoord' and 
 * post-periodic boundary** 'workVector'. Each simulation step repeats this process.  When
 * 'intervalCount' IS zero, 'atomPBIarray' and the atoms' current positions are used as
 * stated below to create the file output.   
 * 
 *   
 * -----'atomPBIarray'----
 * Called in AfterPBC's 'intervalAction', this is an integer array containing a running
 * total of the number of box lengths traversed by the atoms during the simulation.  At 
 * a file write, 'atomPBIarray' is used as a multiplier of the primary box length.  
 * This quantity is added to the atoms' current coordinates to generate an atom's 
 * "absolute displacement" value.  
 * -----------------------
 * 
 *    * and ** : notice the getPriority() methods and their returns
 */


public class MSDCoordWriter implements IAction {

	public MSDCoordWriter(Integrator integrator, Box box, String fileName, int writeInterval) {
		this.box = box;
		iterator = new AtomIteratorLeafAtoms(box);
		afterPBCinstance = new AfterPBC(box, iterator);
		integrator.getEventManager().addListener(new IntegratorListenerAction(afterPBCinstance));
		integrator.getEventManager().addListener(new IntegratorListenerAction(this));

		setWriteInterval(writeInterval);
		double timestep = ((IntegratorMD)integrator).getTimeStep();
		try {
			fileWriter = new FileWriter(fileName, false);
//			fileWriter.write(iterator.size() + "\n");
			fileWriter.write("   write_interval " + writeInterval + "      timestep " + timestep+ "\n");
		} catch (IOException e) {
			System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
		}
	}

	public AfterPBC getAfterPBC() {
		return afterPBCinstance;
	}

	public void setIterator(AtomIterator newIterator) {
        iterator = newIterator;
        afterPBCinstance.setIterator(iterator);
    }

	public void closeFile(){
        try {
            fileWriter.close();
            fileWriter = null;
        }
        catch(IOException e) {
            System.err.println("Cannot close a file, caught IOException: " + e.getMessage());
        }
    }
	
	/**
	 * Equates intervalAction method variable to variable passed in from simulation class
	 * @param writeInterval
	 */
	public void setWriteInterval(int writeInterval){
		this.writeInterval = writeInterval;
		intervalCount = writeInterval;
	}
	
	public void actionPerformed() {
		// our listener may fire after we do; ensure we've unwrapped these coordinates
		afterPBCinstance.unwrap();
		if (--intervalCount == 0){
			Vector boxdim = box.getBoundary().getBoxSize();
			// Gets atomPBIarray from AfterPBC subclass, through the subclass instance
			int [][] atomPBIarray = afterPBCinstance.getAtomPBIarray();

			try {
				iterator.reset();
				int i=0;
				for (IAtom atom = iterator.nextAtom();
                     atom != null; atom = iterator.nextAtom()) {
					Vector atomPosition = atom.getPosition();
					for (int j=0;j < boxdim.getD();j++){
						double actualDistance;
							
						// Total distance traveled between file writes is computed
						actualDistance = atomPBIarray[i][j] * boxdim.getX(j) + atomPosition.getX(j);
						fileWriter.write(""+actualDistance);
						if(j!=boxdim.getD()-1){
							fileWriter.write(" ");
						}
						
					}
					fileWriter.write(" ");
					i++;
				}
				fileWriter.write("\n");
			}
			catch (IOException e) {
	            throw new RuntimeException(e);
	        }
			// Variable is reset after a file write
			intervalCount = writeInterval;
		}
	}

	private AfterPBC afterPBCinstance;
	private Box box;
	private AtomIterator iterator;
	private int writeInterval;
	private int intervalCount;
	private FileWriter fileWriter;

	public static class AfterPBC implements IAction {

		public AfterPBC(Box box, AtomIterator iterator) {
			workVector = box.getSpace().makeVector();
			this.iterator = iterator;
			atomOldCoord = new Vector[box.getLeafList().size()];
			for (int j = 0; j < atomOldCoord.length; j++) {
				atomOldCoord[j] = box.getSpace().makeVector();
			}

			atomPBIarray = new int[box.getLeafList().size()][box.getSpace().D()];
			boxDim = box.getBoundary().getBoxSize();
			interval = 10;
			initAtomOldCoord();
		}

		public void setInterval(int newInterval) {
			interval = newInterval;
			intervalCountdown = interval;
		}
		
		// Method called in main class (see above)
		public int [][] getAtomPBIarray(){
			return atomPBIarray;
		}

		public void setIterator(AtomIterator newIterator) {
            iterator = newIterator;
        }
		
		// Method called in main and sub class (see directly above and above)
		public void initAtomOldCoord() {
			iterator.reset();
			int i=0;
			for (IAtom atom = iterator.nextAtom();
				 atom != null; atom = iterator.nextAtom()) {
				atomOldCoord[i].E(atom.getPosition());
				i++;
			}
			intervalCountdown = interval;
		}

		public void actionPerformed() {
			intervalCountdown--;
			if (intervalCountdown == 0) unwrap();
		}

		public void unwrap() {
			intervalCountdown = interval;
			iterator.reset();
			int i=0;

			// workVector is modified to hold a value of box lengths an atom has traveled
			// atomPBIarray is filled here
			for (IAtom atom = iterator.nextAtom();
				 atom != null; atom = iterator.nextAtom()) {
				workVector.E(atomOldCoord[i]);
				workVector.ME(atom.getPosition());
				workVector.DE(boxDim);

				for (int j = 0; j < boxDim.getD(); j++){

					// Before Math.round, workVector is -/+ 0.9999,1.000,1.0001,0.000
					// Value will truncate when added to atomPBIarray, we must make workVector a whole number
					atomPBIarray[i][j] += Math.round(workVector.getX(j));

				}
				atomOldCoord[i].E(atom.getPosition());
				i++;
			}
		}
		
		private Vector boxDim;
		private int [][] atomPBIarray;
		private Vector workVector;
		private Vector[] atomOldCoord;
		private AtomIterator iterator;
		private int interval, intervalCountdown;
	}
}
