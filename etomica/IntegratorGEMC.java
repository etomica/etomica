package simulate;

import java.util.Random;

public class IntegratorGEMC extends Integrator {
    
    private final Random rand = new Random();
    public double maxRStep, maxVStep;
    public double pressure, betaMu, eBetaMu;
    private int freqDisplace, freqVolume, freqMolecule;
    private int iDisplace, iVolume, iMolecule, iTotal;
    public Phase secondPhase;
    private MCMove atomDisplace = new MCMoveAtom();
    
    public IntegratorGEMC() {
        super();
        freqDisplace = 1;
        freqVolume = freqMolecule = 1;
        maxRStep = 0.1;
        maxVStep = 0.1;
        phaseCountMax = 2;
//        super.phaseCountMax = 2;
        phase = new Phase[phaseCountMax];
        atomDisplace.setAdjustInterval(10000);
        atomDisplace.parentIntegrator = this;
    }
    
  public void registerPhase(Phase p) {
    super.registerPhase(p);
    if(phaseCount > phaseCountMax) {return;}
    if(phaseCount == 2) secondPhase = phase[1];
    System.out.println("phaseCount "+phaseCount);
  }
  
    public void doStep(double dummy) {
        int i = (int)(rand.nextDouble()*iTotal);
        if(i < iDisplace) {
//            if(rand.nextDouble() < 0.5) {trialDisplace(firstPhase);}
//            else {trialDisplace(secondPhase);}
            if(rand.nextDouble() < 0.5) {atomDisplace.doTrial(firstPhase);}
            else {atomDisplace.doTrial(secondPhase);}
        }
        else if(i < iVolume) {
            trialVolume();
        }
        else {
            if(rand.nextDouble() < 0.5) {trialExchange(firstPhase,secondPhase,firstPhase.parentSimulation.firstSpecies);}
            else {trialExchange(secondPhase,firstPhase,firstPhase.parentSimulation.firstSpecies);}
        }
    }
        
    private void trialVolume() {
        double v1Old = firstPhase.volume();
        double v2Old = secondPhase.volume();
        double hOld = firstPhase.potentialEnergy.currentValue() + secondPhase.potentialEnergy.currentValue();
        double vStep = (2.*rand.nextDouble()-1.)*maxVStep;
        double v1New = v1Old + vStep;
        double v2New = v2Old - vStep;
        double v1Scale = v1New/v1Old;
        double v2Scale = v2New/v2Old;
        double r1Scale = Math.pow(v1Scale,1.0/(double)Simulation.D);
        double r2Scale = Math.pow(v2Scale,1.0/(double)Simulation.D);
        firstPhase.inflate(r1Scale);
        secondPhase.inflate(r2Scale);
/*        for(Molecule m=firstPhase.firstMolecule(); m!=null; m=m.nextMolecule()) {  this is done by phase
            m.coordinate.inflate(r1Scale);
        }
        for(Molecule m=secondPhase.firstMolecule(); m!=null; m=m.nextMolecule()) {
            m.coordinate.inflate(r2Scale);
        }*/
        double hNew = firstPhase.potentialEnergy.currentValue() + secondPhase.potentialEnergy.currentValue();
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/temperature+
                       firstPhase.moleculeCount*Math.log(v1Scale) +
                       secondPhase.moleculeCount*Math.log(v2Scale))
                < rand.nextDouble()) 
            {  //reject
              firstPhase.inflate(1.0/r1Scale);
//              for(Molecule m=firstPhase.firstMolecule(); m!=null; m=m.nextMolecule()) {
//                m.coordinate.replace();
//              }
              secondPhase.inflate(1.0/r2Scale);
//              for(Molecule m=secondPhase.firstMolecule(); m!=null; m=m.nextMolecule()) {
//                m.coordinate.replace();
//              }
            }
    }
    
    private void trialExchange(Phase iPhase, Phase dPhase, Species species) {
        
        Species.Agent iSpecies = species.getAgent(iPhase);
        Species.Agent dSpecies = species.getAgent(dPhase);
        double uNew, uOld;
        
        if(dSpecies.nMolecules == 0) {return;}
        
        Molecule m = dSpecies.randomMolecule();
        uOld = dPhase.potentialEnergy.currentValue(m);
        m.displaceTo(iPhase.randomPosition());
        uNew = iPhase.potentialEnergy.insertionValue(m);
        if(uNew == Double.MAX_VALUE) {  //overlap
            m.replace(); 
            return;        
        }        
        double bFactor = dSpecies.nMolecules/dPhase.volume()
                         * iPhase.volume()/(iSpecies.nMolecules+1)
                         * Math.exp(-(uNew-uOld)/temperature);
        if(bFactor > 1.0 || bFactor > rand.nextDouble()) {  //accept
//            dSpecies.deleteMolecule(m);
//            iSpecies.addMolecule(m);
//            dPhase.deleteMolecule(m,dSpecies);
            iPhase.addMolecule(m,iSpecies);  //this handles deletion from dPhase too
        }
        else {              //reject
            m.replace();
        }
    }
           
    public final int getFreqVolume() {return freqVolume;}
    public final void setFreqVolume(int f) {freqVolume = f;}
    public final int getFreqMolecule() {return freqMolecule;}
    public final void setFreqMolecule(int f) {freqMolecule = f;}
    public final int getFreqDisplace() {return freqDisplace;}
    public final void setFreqDisplace(int f) {freqDisplace = f;}
        
    public final double getMaxRStep() {return maxRStep;}
    public final void setMaxRStep(double s) {maxRStep = s;}
    
    public final double getMaxVStep() {return maxVStep;}
    public final void setMaxVStep(double s) {maxVStep = s;}
    
    public void initialize() {
        deployAgents();
        iDisplace = freqDisplace * firstPhase.moleculeCount;
        iVolume = iDisplace + freqVolume;
        iMolecule = iVolume + freqMolecule*firstPhase.moleculeCount;    
        iTotal = iMolecule;
    }
    
    public Integrator.Agent makeAgent(Atom a) {
        return new Agent(a);
    }
    
    public class Agent implements Integrator.Agent {
        public Atom atom;
        public Agent(Atom a) {atom = a;}
    }

}