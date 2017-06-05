/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveSwapConfiguration;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.integrator.mcmove.MCMoveTrialInitiatedEvent;
import etomica.space.Space;
import etomica.units.Null;
import etomica.util.IEvent;
import etomica.util.IListener;

/**
 * Parallel-tempering integrator.  Oversees other integrators that are defined to perform
 * MC trials (or perhaps molecular dynamics) in different boxes.  These integrators
 * are identified (added) to this integrator, and the simulation runs on this
 * integrator's thread.  When this integrator does a step, it passes on the
 * instruction to each added integrator, causing all to do a step on each's box.
 * Occasionally, this integrator will instead attempt to swap the configurations
 * of two of the boxes.  Acceptance of this move depends on parameters (usually 
 * temperature) of the integrators for the boxes, as well as the configurations
 * in each box.  The swap trial is performed by a MCMove class designed for
 * this purpose.  Such a class is made by a MCMoveSwapConfigurationFactory class
 * that is identified to this integrator (a default is selected if not specified).
 * Every time an integrator is added to this one, a MCMoveSwap class is made (by this
 * integrator using the factory) to manage swap trials between the new integrator's
 * box and that of the one most recently added.
 * 
 * @author David Kofke
 */
public class IntegratorPT extends IntegratorManagerMC {
    
    private final Space space;

    public IntegratorPT(IRandom random, Space _space) {
        this(random, MCMoveSwapConfiguration.FACTORY, _space);
    }
    
    public IntegratorPT(IRandom random, MCMoveSwapFactory swapFactory, Space _space) {
        super(random);
        this.space = _space;
        setGlobalMoveInterval(100);
        mcMoveSwapFactory = swapFactory;
    }

    /**
     * Adds the given integrator to those managed by this integrator, and
     * includes integrator's box to the set among which configurations are
     * swapped.  Box of new integrator will be subject to swapping with
     * the most recently added integrator/box (and the next one, if another
     * is added subsequently).
     */
	public void addIntegrator(IntegratorBox integrator){
	    super.addIntegrator(integrator);
        
		if (nIntegrators > 1) {
            MCMove newMCMove = mcMoveSwapFactory.makeMCMoveSwap((IntegratorBox)integrators[nIntegrators-2], 
                    (IntegratorBox)integrators[nIntegrators-1], space);
            moveManager.addMCMove(newMCMove);
		}
	}
    
    private static final long serialVersionUID = 1L;
	private final MCMoveSwapFactory mcMoveSwapFactory;

	
	/**
	 * Interface for a class that can make a MCMove that will swap
	 * the configurations of two boxes.  Different MCMove classes
	 * would do this differently, depending on ensemble of simulation
	 * and other factors.
	 */
	public interface MCMoveSwapFactory {
	    /**
	     * @param integrator1 integrator for one of the boxes being swapped
	     * @param integrator2 integrator for the other box
	     */
	    public MCMove makeMCMoveSwap(IntegratorBox integrator1, IntegratorBox integrator2, Space _space);
	}

    /**
     * Interface for a move that swaps two boxes.  Enables access to
     * the swapped boxes.
     */
    public interface MCMoveSwap {
        public Box[] swappedBoxes();
    }


	/**
     * Meter that tracks the swapping of the boxes in parallel-tempering
     * simulation.  Designed for input to a DisplayPlot to provide a graphical
     * record of how the boxes swap configurations.
     */
    public static class BoxTracker implements IEtomicaDataSource, IListener, java.io.Serializable {
        
        public BoxTracker() {
            data = new DataDoubleArray(0);
            dataInfo = new DataInfoDoubleArray("Box Tracker", Null.DIMENSION, new int[]{0});
            tag = new DataTag();
            dataInfo.addTag(tag);
        }
        
        public IEtomicaDataInfo getDataInfo() {
            return dataInfo;
        }
        
        public DataTag getTag() {
            return tag;
        }
        
        /**
         * Method called when two boxes are successfully exchanged.
         */
        public void actionPerformed(IEvent evt) {
            if(evt instanceof MCMoveTrialInitiatedEvent || !((MCMoveTrialCompletedEvent)evt).isAccepted()) return;
            if(!(((MCMoveEvent)evt).getMCMove() instanceof MCMoveSwap)) return;
            Box[] boxes = ((MCMoveSwap)((MCMoveEvent)evt).getMCMove()).swappedBoxes();
            int i0 = boxes[0].getIndex()-1;
            int i1 = boxes[1].getIndex()-1;
            int temp = track[i0];
            track[i0] = track[i1];
            track[i1] = temp;
        }
        
        /**
         * Specifies the number of boxes that are tracked.
         */
        public void setNumBoxes(int numBoxes) {
            track = new int[numBoxes];
            data = new DataDoubleArray(new int[] {numBoxes});
            dataInfo = new DataInfoDoubleArray("Box Tracker", Null.DIMENSION, new int[]{numBoxes});
            dtrack = data.getData();
            for(int i=0; i<numBoxes; i++) {
                track[i] = i;
            }
        }
        
        /**
         * Returns array y such that y[i] is the current
         * box of the configuration that began in box i.
         */
        public IData getData() {
            for (int i=0; i<track.length; i++) {
                dtrack[i] = track[i];
            }
            return data;
        }
        
        public int getDataLength() {
            return track.length;
        }
        
        private static final long serialVersionUID = 1L;
        private int[] track;
        private double[] dtrack;
        private DataDoubleArray data;
        private DataInfoDoubleArray dataInfo;
        private final DataTag tag;
    }
    
}

