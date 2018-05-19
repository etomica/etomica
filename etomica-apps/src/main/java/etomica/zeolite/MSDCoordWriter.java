/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.zeolite;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.IAction;
import etomica.action.activity.ControllerEvent;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.space.Vector;
import etomica.atom.iterator.AtomIteratorBoxDependent;
import etomica.integrator.IntegratorListenerAction;
import etomica.space.Space;
import etomica.util.IEvent;
import etomica.util.IListener;

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


public class MSDCoordWriter implements IAction, IListener {
	
	public MSDCoordWriter(Space _space, String fileName) throws RuntimeException {
	    throw new RuntimeException("MSDCoordWriter is not usable due to the removal of " +
	                               "action priorities.");
/*
		// Creates an instance of subclass AfterPBC
		iterator = new AtomIteratorLeafAtoms();
        afterPBCinstance = new AfterPBC(_space,iterator);
		this.fileName = fileName;
		setWriteInterval(1);
*/
	}
	
	public void setBox(Box newBox){
		
		box = newBox;
		iterator.setBox(box);
		afterPBCinstance.setBox(box);
	}
    
    public void setIterator(AtomIteratorBoxDependent newIterator) {
        iterator = newIterator;
        afterPBCinstance.setIterator(iterator);
    }
	
	public void setIntegrator(Integrator integrator){
		integrator.getEventManager().addListener(new IntegratorListenerAction(this));
//        integrator.setIntervalActionPriority(this, 50);
		integrator.getEventManager().addListener(new IntegratorListenerAction(afterPBCinstance));
//        integrator.setIntervalActionPriority(afterPBCinstance, 200);
	}
	
	// Methods involved with file creation/closing
	public void openFile(){
		try { 
			fileWriter = new FileWriter(fileName, false);
			fileWriter.write(iterator.size()+"\n");
		}
	    catch(IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
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
		afterPBCinstance.updateAtomOldCoord();
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
							fileWriter.write("\t");
						}
						
					}
					fileWriter.write("\n");
					i++;
				}
			}
			catch (IOException e) {
	            throw new RuntimeException(e);
	        }
			// Variable is reset after a file write
			intervalCount = writeInterval;
		}
	}
	
	// *
	public int getPriority() {
		return 50;
	}

    public void actionPerformed(IEvent evt) {
        if (fileWriter != null &&
            (((ControllerEvent)evt).getType() == ControllerEvent.Type.NO_MORE_ACTIONS ||
             ((ControllerEvent)evt).getType() == ControllerEvent.Type.HALTED)) {
            closeFile();
        }
    }

	private AfterPBC afterPBCinstance;
	private Box box;
	private AtomIteratorBoxDependent iterator;
	private int writeInterval;
	private int intervalCount;
	private String fileName;
	private FileWriter fileWriter;
	
	/*
	 * -------------------------SUBCLASS AfterPBC----------------------------------
	 */
	
	private static class AfterPBC implements IAction {
		
		public AfterPBC(Space _space, AtomIteratorBoxDependent iterator){
			workVector = _space.makeVector();
			this.iterator = iterator;
			this.space = _space;
		}
		
		// Method called in main class (see above)
		public int [][] getAtomPBIarray(){
			return atomPBIarray;
		}
				
		public void setBox(Box box){
			atomOldCoord = new Vector[box.getLeafList().size()];
			for (int j=0; j < atomOldCoord.length; j++){
				atomOldCoord[j] = space.makeVector();
			}
						
			atomPBIarray = new int[box.getLeafList().size()][space.D()];
			iterator.setBox(box);
			boxDim = box.getBoundary().getBoxSize();
			updateAtomOldCoord();
		}
        
        public void setIterator(AtomIteratorBoxDependent newIterator) {
            iterator = newIterator;
        }
		
		// Method called in main and sub class (see directly above and above)
		public void updateAtomOldCoord(){
			iterator.reset();
			int i=0;
            for (IAtom atom = iterator.nextAtom();
                 atom != null; atom = iterator.nextAtom()) {
				atomOldCoord[i].E(atom.getPosition());
				i++;
			}
		}
		
		public void actionPerformed() {
			iterator.reset();
			int i=0;
			
			// workVector is modified to hold a value of box lengths an atom has traveled
			// atomPBIarray is filled here
            for (IAtom atom = iterator.nextAtom();
                 atom != null; atom = iterator.nextAtom()) {
				workVector.E(atomOldCoord[i]);
				workVector.ME(atom.getPosition()); 
				workVector.DE(boxDim);
				
				for (int j=0;j < boxDim.getD();j++){
					
					// Before Math.round, workVector is -/+ 0.9999,1.000,1.0001,0.000
					// Value will truncate when added to atomPBIarray, we must make workVector a whole number
					atomPBIarray[i][j] += Math.round(workVector.getX(j));
					
				}
				i++;
			}
		}
		
		private Vector boxDim;
		private int [][] atomPBIarray;
		private Vector workVector;
		private Vector[] atomOldCoord;
		private AtomIteratorBoxDependent iterator;
		private final Space space;
	}
	
}
