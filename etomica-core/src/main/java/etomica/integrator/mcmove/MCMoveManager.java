/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import java.util.ArrayList;
import java.util.List;

import etomica.box.Box;
import etomica.meta.annotations.IgnoreProperty;
import etomica.util.random.IRandom;

public class MCMoveManager {

    public MCMoveManager(IRandom random) {
        this.random = random;
        mcMoveList = new ArrayList<MCMoveLinker>();
    }

    /**
     * Sets moves in given array to be integrator's set of moves, deleting any
     * existing moves.
     */
    public void setMCMoves(MCMove[] moves) {
        mcMoveList.clear();
        for (int i = 0; i < moves.length; i++) {
            addMCMove(moves[i]);
        }
    }

    /**
     * Constructs and returns array of all moves added to the integrator.
     */
    public List<MCMove> getMCMoves() {
        List<MCMove> l = new ArrayList<MCMove>();
        for (int i=0; i<mcMoveList.size(); i++) {
            l.add(mcMoveList.get(i).move);
        }
        return l;
    }

    /**
     * Adds the given MCMove to the set of moves performed by the integrator and
     * recalculates move frequencies.  If the MCMoveManager has been given a 
     * Box, the MCMove added here must be of type MCMoveBox.
     * 
     * @throws ClassCastException if this MCMoveManager has a Box and the
     * given MCMove is not an MCMoveBox
     */
    public void addMCMove(MCMove move) {
        //make sure move wasn't added already
        for (int i=0; i<mcMoveList.size(); i++) {
            if (mcMoveList.get(i).move == move) throw new RuntimeException("move already added");
        }
        mcMoveList.add(new MCMoveLinker(move));

        MCMoveTracker tracker = move.getTracker();
        if (tracker instanceof MCMoveStepTracker && !isEquilibrating) {
            ((MCMoveStepTracker)tracker).setTunable(isEquilibrating);
        }
        if (box != null) {
            ((MCMoveBox)move).setBox(box);
        }
        recomputeMoveFrequencies();
    }

    /**
     * Removes the given MCMove from the set of moves performed by the integrator and
     * recalculates move frequencies.  Returns false if the move was not used by
     * the integrator.
     */
    public boolean removeMCMove(MCMove move) {
        //make sure move wasn't added already
        int pos = -1;
        for (int i=0; i<mcMoveList.size(); i++) {
            if (mcMoveList.get(i).move == move) {
                pos = i;
                break;
            }
        }
        if (pos == -1) return false;
        mcMoveList.remove(pos);
        recomputeMoveFrequencies();
        return true;
    }

    /**
     * Invokes superclass method and informs all MCMoves about the new box.
     * The moves are assumed to all be of type MCMoveBox.
     * 
     * @throws ClassCastException if any move is not an MCMoveBox
     */
    public void setBox(Box p) {
        box = p;
        for (int i=0; i<mcMoveList.size(); i++) {
            ((MCMoveBox)mcMoveList.get(i).move).setBox(box);
        }
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    /**
     * Selects a MCMove instance from among those added to the integrator, with
     * probability in proportion to the frequency value assigned to the move.
     */
    public MCMove selectMove() {
//        System.out.println("select move");
        if (mcMoveList.size() == 0) return null;
        int i = random.nextInt(frequencyTotal);
        int pos = 0;
        while ((i -= mcMoveList.get(pos).fullFrequency) >= 0) {
            pos++;
        }
        selectedLink = mcMoveList.get(pos);
        return selectedLink.move;
    }

    @IgnoreProperty
    public MCMove getSelectedMove() {
        return selectedLink.move;
    }

    /**
     * Recomputes all the move frequencies.
     */
    public void recomputeMoveFrequencies() {
        frequencyTotal = 0;
        for (int i=0; i<mcMoveList.size(); i++) {
            MCMoveLinker link = mcMoveList.get(i);
            link.resetFullFrequency();
            frequencyTotal += link.fullFrequency;
        }
    }
    
    public void setEquilibrating(boolean equilibrating) {
        isEquilibrating = equilibrating;
        for (int i=0; i<mcMoveList.size(); i++) {
            MCMoveTracker tracker = mcMoveList.get(i).move.getTracker();
            if (tracker instanceof MCMoveStepTracker) {
                ((MCMoveStepTracker)tracker).setTunable(isEquilibrating);
            }
        }
    }

    /**
     * Returns the trial frequency set for the given move, over and above the
     * move's nominal frequency.  The frequency is 1.0 by default.
     * 
     * This method should not be called with a move not contained by this
     * manager.
     */
    public double getFrequency(MCMove move) {
        for (int i=0; i<mcMoveList.size(); i++) {
            MCMoveLinker link = mcMoveList.get(i);
            if (link.move == move) {
                return link.frequency;
            }
        }
        throw new IllegalArgumentException("I don't have "+move);
    }

    /**
     * Sets the trial frequency for the given move to the given frequency.  The
     * frequency is used in conjunction with the move's own nominal frequency.
     * 
     * This method should not be called with a move not contained by this
     * manager.
     */
    public void setFrequency(MCMove move, double newFrequency) {
        for (int i=0; i<mcMoveList.size(); i++) {
            MCMoveLinker link = mcMoveList.get(i);
            if (link.move == move) {
                link.frequency = newFrequency;
                recomputeMoveFrequencies();
                return;
            }
        }
        throw new IllegalArgumentException("I don't have "+move);
    }

    private Box box;
    protected final List<MCMoveLinker> mcMoveList;
    private MCMoveLinker selectedLink;
    private int frequencyTotal;
    private boolean isEquilibrating = true;
    private final IRandom random;

    /**
     * Linker used to construct linked-list of MCMove instances
     */
    private static class MCMoveLinker {
        protected int fullFrequency;
        protected double frequency;
        protected final MCMove move;

        MCMoveLinker(MCMove move) {
            this.move = move;
            frequency = 1.0;
        }

        /**
         * Updates the full frequency based on the current value of the
         * frequency, the status of the perParticleFrequency flag, and the
         * current number of molecules in the boxes affected by the move.
         */
        void resetFullFrequency() {
            fullFrequency = (int)Math.round(frequency * move.getNominalFrequency());
            if ((move instanceof MCMoveBox) && ((MCMoveBox)move).getBox() != null
                    && ((MCMoveBox)move).isNominallyPerParticleFrequency() ) {
                int nMolecules = ((MCMoveBox)move).getBox().getMoleculeList().size();
                if (nMolecules > 0) fullFrequency *= nMolecules;
            }
        }

    }
}
