package simulate;

import java.util.Random;

public class IntegratorGEMC extends Integrator {
    
    private final Random rand = new Random();
    public double maxRStep, maxVStep;
    public double pressure, betaMu, eBetaMu;
    private int freqDisplace, freqVolume, freqMolecule;
    private int iDisplace, iVolume, iMolecule, iTotal;
    public PhaseSpace secondPhaseSpace;
    private MCMove atomDisplace = new MCMoveAtom();
    
    public IntegratorGEMC() {
        super();
        freqDisplace = 1;
        freqVolume = freqMolecule = 1;
        maxRStep = 0.1;
        maxVStep = 0.1;
        phaseSpaceCountMax = 2;
//        super.phaseSpaceCountMax = 2;
        phaseSpace = new PhaseSpace[phaseSpaceCountMax];
        atomDisplace.setAdjustInterval(10000);
        atomDisplace.parentIntegrator = this;
    }
    
  public void registerPhaseSpace(PhaseSpace p) {
    super.registerPhaseSpace(p);
    if(phaseSpaceCount > phaseSpaceCountMax) {return;}
    if(phaseSpaceCount == 2) secondPhaseSpace = phaseSpace[1];
    System.out.println("phaseSpaceCount "+phaseSpaceCount);
  }
  
    public void doStep(double dummy) {
        int i = (int)(rand.nextDouble()*iTotal);
        if(i < iDisplace) {
//            if(rand.nextDouble() < 0.5) {trialDisplace(firstPhase);}
//            else {trialDisplace(secondPhase);}
            if(rand.nextDouble() < 0.5) {atomDisplace.doTrial(firstPhaseSpace);}
            else {atomDisplace.doTrial(secondPhaseSpace);}
        }
        else if(i < iVolume) {
            trialVolume();
        }
        else {
            if(rand.nextDouble() < 0.5) {trialExchange(firstPhaseSpace.firstSpecies(),secondPhaseSpace.firstSpecies());}
            else {trialExchange(secondPhaseSpace.firstSpecies(),firstPhaseSpace.firstSpecies());}
        }
    }
        
/*    private void trialDisplace(Phase p) {
        if(p.nMoleculeTotal == 0) {return;}
        double uOld, uNew;
        int i = (int)(rand.nextDouble()*p.nAtomTotal);
        AtomC a = (AtomC)p.firstAtom();
        for(int j=i; --j>=0; ) {a = a.getNextAtomC();}  //get ith atom in list
        uOld = p.potentialEnergy.currentValue(a);
        Space.randomVector(dr, maxRStep, rand);
        a.displace(dr);
        uNew = p.potentialEnergy.currentValue(a);
        if(uNew < uOld) {return;}   //accept
        if(uNew >= Double.MAX_VALUE || 
           Math.exp(-(uNew-uOld)/temperature) < rand.nextDouble()) { //reject
             a.replace();
             return;
        }           
    }
    */
    private void trialVolume() {
        double v1Old = firstPhaseSpace.volume();
        double v2Old = secondPhaseSpace.volume();
        double hOld = firstPhaseSpace.potentialEnergy.currentValue() + secondPhaseSpace.potentialEnergy.currentValue();
        double vStep = (2.*rand.nextDouble()-1.)*maxVStep;
        double v1New = v1Old + vStep;
        double v2New = v2Old - vStep;
        double v1Scale = v1New/v1Old;
        double v2Scale = v2New/v2Old;
        double r1Scale = Math.pow(v1Scale,1.0/(double)Simulation.D);
        double r2Scale = Math.pow(v2Scale,1.0/(double)Simulation.D);
        firstPhaseSpace.inflate(r1Scale);
        secondPhaseSpace.inflate(r2Scale);
        for(Molecule m=firstPhaseSpace.firstMolecule(); m!=null; m=m.nextMolecule()) {
            m.coordinate.inflate(r1Scale);
        }
        for(Molecule m=secondPhaseSpace.firstMolecule(); m!=null; m=m.nextMolecule()) {
            m.coordinate.inflate(r2Scale);
        }
        double hNew = firstPhaseSpace.potentialEnergy.currentValue() + secondPhaseSpace.potentialEnergy.currentValue();
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/temperature+
                       firstPhaseSpace.moleculeCount*Math.log(v1Scale) +
                       secondPhaseSpace.moleculeCount*Math.log(v2Scale))
                < rand.nextDouble()) 
            {  //reject
              firstPhaseSpace.inflate(1.0/r1Scale);
              for(Molecule m=firstPhaseSpace.firstMolecule(); m!=null; m=m.nextMolecule()) {
                m.coordinate.replace();
              }
              secondPhaseSpace.inflate(1.0/r2Scale);
              for(Molecule m=secondPhaseSpace.firstMolecule(); m!=null; m=m.nextMolecule()) {
                m.coordinate.replace();
              }
            }
    }
    
    private void trialExchange(Species iSpecies, Species dSpecies) {
        double uNew, uOld;
        
        if(dSpecies.nMolecules == 0) {return;}
        
        Molecule m = dSpecies.randomMolecule();
        uOld = dSpecies.parentPhaseSpace.potentialEnergy.currentValue(m);
        m.coordinate.displaceToRandom(iSpecies.parentPhaseSpace.dimensions());
//        m.parentSpecies = iSpecies;
        uNew = iSpecies.parentPhaseSpace.potentialEnergy.insertionValue(m);
        if(uNew == Double.MAX_VALUE) {
            m.replace(); 
//            m.parentSpecies = dSpecies; 
            return;        //overlap
        }        
        double bFactor = dSpecies.nMolecules/dSpecies.parentPhaseSpace.volume()
                         * iSpecies.parentPhaseSpace.volume()/(iSpecies.nMolecules+1)
                         * Math.exp(-(uNew-uOld)/temperature);
        if(bFactor > 1.0 || bFactor > rand.nextDouble()) {  //accept
//            m.parentSpecies = dSpecies; 
            dSpecies.deleteMolecule(m);
            iSpecies.addMolecule(m);
        }
        else {              //reject
            m.replace();
//            m.parentSpecies = dSpecies;
            return;
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
        iDisplace = freqDisplace * firstPhaseSpace.moleculeCount;
        iVolume = iDisplace + freqVolume;
        iMolecule = iVolume + freqMolecule*firstPhaseSpace.moleculeCount;    
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