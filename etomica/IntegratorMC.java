package simulate;

import java.util.Random;

public class IntegratorMC extends Integrator {
    
    private final Random rand = new Random();
    private MCMove firstMove, lastMove;
    private int frequencyTotal;
    private transient MCMove trialMove;
    
    public IntegratorMC() {
        super();
    }
    
    public void add(MCMove move) {
        if(firstMove == null) {firstMove = move;}
        else {lastMove.setNextMove(move);}
        lastMove = move;
        move.parentIntegrator = this;
    }
                
    public void doStep(double dummy) {
        int i = (int)(rand.nextDouble()*frequencyTotal);
        trialMove = firstMove;
        while((i-=trialMove.getFrequency()) >= 0) {
            trialMove = trialMove.getNextMove();
        }
        trialMove.doTrial(firstPhase);
    }
    
    public void initialize() {
        frequencyTotal = 0;
        for(MCMove m=firstMove; m!=null; m=m.getNextMove()) {
            m.resetFrequency(firstPhase);
            frequencyTotal += m.getFrequency();
        }
    }
}