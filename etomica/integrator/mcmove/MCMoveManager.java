package etomica.integrator.mcmove;

import java.io.Serializable;

import etomica.phase.Phase;
import etomica.simulation.Simulation;

public class MCMoveManager implements Serializable {

    public MCMoveManager() {
        super();
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
     * recalculates move frequencies.
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
        if (phase != null && move instanceof MCMovePhase) {
            ((MCMovePhase)move).setPhase(phase);
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
                link.nextLink = link.nextLink.nextLink;
                moveCount--;
                recomputeMoveFrequencies();
                return true;
            }
        }
        return false;
    }

    /**
     * Invokes superclass method and informs all MCMoves about the new phase.
     * Moves are not notified if they have a number of phases different from
     * the number of phases handled by the integrator.
     */
    public void setPhase(Phase p) {
        phase = p;
        for (MCMoveLinker link = firstMoveLink; link != null; link = link.nextLink) {
            if (link.move instanceof MCMovePhase) {
                ((MCMovePhase)link.move).setPhase(phase);
            }
        }
    }

    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }

    /**
     * Selects a MCMove instance from among those added to the integrator, with
     * probability in proportion to the frequency value assigned to the move.
     */
    public MCMove selectMove() {
        if (firstMoveLink == null)
            return null;
        int i = Simulation.random.nextInt(frequencyTotal);
        MCMoveLinker link = firstMoveLink;
        while ((i -= link.fullFrequency) >= 0) {
            link = link.nextLink;
        }
        link.selectionCount++;
        return link.move;
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

    private Phase phase;
    private MCMoveLinker firstMoveLink, lastMoveLink;
    private int frequencyTotal;
    private int moveCount;
    private boolean isEquilibrating;

    /**
     * Linker used to construct linked-list of MCMove instances
     */
    private static class MCMoveLinker implements java.io.Serializable {
        int frequency, fullFrequency;
        final MCMove move;
        boolean perParticleFrequency;
        int selectionCount;
        MCMoveLinker nextLink;

        MCMoveLinker(MCMove move) {
            this.move = move;
            frequency = move.getNominalFrequency();
            perParticleFrequency = (move instanceof MCMovePhase) && ((MCMovePhase)move).isNominallyPerParticleFrequency();
        }

        /**
         * Updates the full frequency based on the current value of the
         * frequency, the status of the perParticleFrequency flag, and the
         * current number of molecules in the phases affected by the move.
         */
        void resetFullFrequency() {
            fullFrequency = frequency;
            if (perParticleFrequency && ((MCMovePhase)move).getPhase() != null) {
                fullFrequency *= ((MCMovePhase)move).getPhase().moleculeCount();
            }
        }

    }
}
