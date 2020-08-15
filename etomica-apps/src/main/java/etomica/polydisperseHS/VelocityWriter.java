/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.polydisperseHS;

import etomica.action.IAction;
import etomica.action.activity.ControllerEvent;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
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

public class VelocityWriter implements IAction {

	public VelocityWriter(Integrator integrator, Box box, String fileName, int writeInterval) {
		this.box = box;
		iterator = new AtomIteratorLeafAtoms(box);
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


	public void setIterator(AtomIterator newIterator) {
        iterator = newIterator;
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
		if (--intervalCount == 0){
			Vector boxdim = box.getBoundary().getBoxSize();
			// Gets atomPBIarray from AfterPBC subclass, through the subclass instance

			try {
				iterator.reset();
				int i=0;
				for (IAtom atom = iterator.nextAtom();
                     atom != null; atom = iterator.nextAtom()) {
					Vector atomVelocity = ((IAtomKinetic)atom).getVelocity();
					for (int j=0;j < boxdim.getD();j++){
						fileWriter.write(""+atomVelocity.getX(j));
						if(j!=boxdim.getD()-1){
							fileWriter.write(" ");
						}
						
					}
					fileWriter.write(" " );
					i++;
				}
				fileWriter.write("\n" );
			}
			catch (IOException e) {
	            throw new RuntimeException(e);
	        }
			// Variable is reset after a file write
			intervalCount = writeInterval;
		}
	}


	private Box box;
	private AtomIterator iterator;
	private int writeInterval;
	private int intervalCount;
	private FileWriter fileWriter;

}