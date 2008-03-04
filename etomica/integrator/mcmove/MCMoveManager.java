package etomica.integrator.mcmove;

import java.io.Serializable;

import etomica.api.IBox;
import etomica.api.IRandom;

public class MCMoveManager implements Serializable {

    public MCMoveManager(IRandom random) {
        this.random = random;
    }

    /**
     * Sets moves in given array to be integrator's set of moves, deleting any
     * existing moves.
     */
    public void setMCMoves(MCMove[] moves) {
        firstMoveLink = null;
        moveCount = 0;
        for (int i = 0; i < moves.length; i++) {
            addMCMove(moves[i]);
        }
    }

    /**
     * Constructs and returns array of all moves added to the integrator.
     */
    public MCMove[] getMCMoves() {
        MCMove[] moves = new MCMove[moveCount];
        int i = 0;
        for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
            moves[i++] = link.move;
        }
        return moves;
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
        for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
            if (move == link.move)
                return;
        }
        if (firstMoveLink == null) {
            firstMoveLink = new MCMoveLinker(move);
            lastMoveLink = firstMoveLink;
        } else {
            lastMoveLink.nextLink = new MCMoveLinker(move);
            lastMoveLink = lastMoveLink.nextLink;
        }
        if (box != null) {
            ((MCMoveBox)move).setBox(box);
        }
        moveCount++;
        recomputeMoveFrequencies();
    }

    /**
     * Removes the given MCMove from the set of moves performed by the integrator and
     * recalculates move frequencies.  Returns false if the move was not used by
     * the integrator.
     */
    public boolean removeMCMove(MCMove move) {
        //make sure move wasn't added already
        if (move == firstMoveLink.move) {
            firstMoveLink = firstMoveLink.nextLink;
            moveCount--;
            recomputeMoveFrequencies();
            return true;
        }
        for (MCMoveLinker link = firstMoveLink; link.nextLink != null; link = link.nextLink) {
            if (move == link.nextLink.move) {
                if (link.nextLink == lastMoveLink) {
                    lastMoveLink = link;
                }
                link.nextLink = link.nextLink.nextLink;
                moveCount--;
                recomputeMoveFrequencies();
                return true;
            }
        }
        return false;
    }

    /**
     * Invokes superclass method and informs all MCMoves about the new box.
     * The moves are assumed to all be of type MCMoveBox.
     * 
     * @throws ClassCastException if any move is not an MCMoveBox
     */
    public void setBox(IBox p) {
        box = p;
        for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
            ((MCMoveBox)link.move).setBox(box);
        }
    }

    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }

    /**
     * Selects a MCMove instance from among those added to the integrator, with
     * probability in proportion to the frequency value assigned to the move.
     */
    public MCMove selectMove() {
        if (firstMoveLink == null || frequencyTotal == 0) {
            selectedLink = null;
            return null;
        }
        int i = random.nextInt(frequencyTotal);
        selectedLink = firstMoveLink;
        while ((i -= selectedLink.fullFrequency) >= 0) {
            selectedLink = selectedLink.nextLink;
        }
        selectedLink.selectionCount++;
        return selectedLink.move;
    }
    
    public MCMove getSelectedMove() {
        return selectedLink.move;
    }
    /**
     * Recomputes all the move frequencies.
     */
    public void recomputeMoveFrequencies() {
        frequencyTotal = 0;
        for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
            link.resetFullFrequency();
            frequencyTotal += link.fullFrequency;
        }
    }
    
    public void setEquilibrating(boolean equilibrating) {
        isEquilibrating = equilibrating;
        for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
            MCMoveTracker tracker = link.move.getTracker();
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
        for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
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
        for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
            if (link.move == move) {
                link.frequency = newFrequency;
                recomputeMoveFrequencies();
                return;
            }
        }
        throw new IllegalArgumentException("I don't have "+move);
    }

    private static final long serialVersionUID = 1L;
    private IBox box;
    private MCMoveLinker firstMoveLink, lastMoveLink;
    private MCMoveLinker selectedLink;
    private int frequencyTotal;
    private int moveCount;
    private boolean isEquilibrating;
    private final IRandom random;

    /**
     * Linker used to construct linked-list of MCMove instances
     */
    private static class MCMoveLinker implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        int fullFrequency;
        double frequency;
        final MCMove move;
        int selectionCount;
        MCMoveLinker nextLink;

        MCMoveLinker(MCMove move) {
            this.move = move;
            frequency = 1.0;
        }

        /**
         * Updates the full frequency based on the current value of the
         * frequency, the status of the perParticleFrequency flag, and the
         * current number of molecules in the boxs affected by the move.
         */
        void resetFullFrequency() {
            fullFrequency = (int)(frequency * move.getNominalFrequency());
            if ((move instanceof MCMoveBox) && ((MCMoveBox)move).getBox() != null
                    && ((MCMoveBox)move).isNominallyPerParticleFrequency() ) {
                fullFrequency *= ((MCMoveBox)move).getBox().moleculeCount();
            }
        }

    }
}
