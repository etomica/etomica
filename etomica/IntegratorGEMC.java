package simulate;

import java.util.Random;

public class IntegratorGEMC extends Integrator {
    
    private transient final double[] dr = new double[Space.D];
    private final Random rand = new Random();
    public double maxRStep, maxVStep;
    public double pressure, betaMu, eBetaMu;
    private int freqDisplace, freqVolume, freqMolecule;
    private int iDisplace, iVolume, iMolecule, iTotal;
    public Phase secondPhase;
    
    public IntegratorGEMC() {
        super();
        freqDisplace = 1;
        freqVolume = freqMolecule = 1;
        maxRStep = 0.1;
        maxVStep = 0.1;
        nPhasesMax = 2;
        super.nPhasesMax = 2;
        phase = new Phase[nPhasesMax];
    }
    
  public void registerPhase(Phase p) {
    super.registerPhase(p);
    if(nPhases > nPhasesMax) {return;}
    if(nPhases == 2) secondPhase = phase[1];
    System.out.println("nPhases "+nPhases);
  }
  
    public void doStep(double dummy) {
        int i = (int)(rand.nextDouble()*iTotal);
        if(i < iDisplace) {
            if(rand.nextDouble() < 0.5) {trialDisplace(firstPhase);}
            else {trialDisplace(secondPhase);}
        }
        else if(i < iVolume) {
            trialVolume();
        }
        else {
            if(rand.nextDouble() < 0.5) {trialExchange(firstPhase.firstSpecies(),secondPhase.firstSpecies());}
            else {trialExchange(secondPhase.firstSpecies(),firstPhase.firstSpecies());}
        }
    }
        
    private void trialDisplace(Phase p) {
        if(p.nMoleculeTotal == 0) {return;}
        double uOld, uNew;
        int i = (int)(rand.nextDouble()*p.nAtomTotal);
        AtomC a = (AtomC)p.firstAtom();
        for(int j=i; --j>=0; ) {a = a.getNextAtomC();}  //get ith atom in list
        uOld = a.potentialEnergy();
        Space.randomVector(dr, maxRStep, rand);
        a.displace(dr);
        uNew = a.potentialEnergy();
        if(uNew < uOld) {return;}   //accept
        if(uNew >= Double.MAX_VALUE || 
           Math.exp(-(uNew-uOld)/temperature) < rand.nextDouble()) { //reject
             a.replace();
             return;
        }           
    }
    
    private void trialVolume() {
        double v1Old = firstPhase.space.volume;
        double v2Old = secondPhase.space.volume;
        double hOld = energy(firstPhase) + energy(secondPhase);
        double vStep = (2.*rand.nextDouble()-1.)*maxVStep;
        double v1New = v1Old + vStep;
        double v2New = v2Old - vStep;
        double v1Scale = v1New/v1Old;
        double v2Scale = v2New/v2Old;
        double r1Scale = Math.pow(v1Scale,1.0/(double)Space.D);
        double r2Scale = Math.pow(v2Scale,1.0/(double)Space.D);
        firstPhase.space.inflate(r1Scale);
        secondPhase.space.inflate(r2Scale);
        for(Molecule m=firstPhase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
            m.inflate(r1Scale);
        }
        for(Molecule m=secondPhase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
            m.inflate(r2Scale);
        }
        double hNew = energy(firstPhase) + energy(secondPhase);
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/temperature+
                       firstPhase.nMoleculeTotal*Math.log(v1Scale) +
                       secondPhase.nMoleculeTotal*Math.log(v2Scale))
                < rand.nextDouble()) 
            {  //reject
              firstPhase.space.inflate(1.0/r1Scale);
              for(Molecule m=firstPhase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
                m.replace();
              }
              secondPhase.space.inflate(1.0/r2Scale);
              for(Molecule m=secondPhase.firstMolecule(); m!=null; m=m.getNextMolecule()) {
                m.replace();
              }
            }
    }
    
    private void trialExchange(Species iSpecies, Species dSpecies) {
        double uNew, uOld;
        
        if(dSpecies.nMolecules == 0) {return;}
        
        Molecule m = dSpecies.randomMolecule();
        uOld = m.potentialEnergy();
        Space.uEv1(dr,m.COM());
        Space.uTEa1(dr,-1.0);
        m.displace(dr);  //zero center-of-mass, saving old coordinates
        iSpecies.parentPhase.space.randomVector(dr, rand);  //random point in insertion volume
        m.translate(dr);
        m.parentSpecies = iSpecies;
        uNew = m.potentialEnergy();
        if(uNew == Double.MAX_VALUE) {
            m.replace(); 
            m.parentSpecies = dSpecies; 
            return;        //overlap
        }        
        double bFactor = dSpecies.nMolecules/dSpecies.parentPhase.space.volume
                         * iSpecies.parentPhase.space.volume/(iSpecies.nMolecules+1)
                         * Math.exp(-(uNew-uOld)/temperature);
        if(bFactor > 1.0 || bFactor > rand.nextDouble()) {  //accept
            m.parentSpecies = dSpecies; 
            dSpecies.deleteMolecule(m);
            iSpecies.addMolecule(m);
        }
        else {              //reject
            m.replace();
            m.parentSpecies = dSpecies;
            return;
        }
    }
           
    //Computes total energy of phase
    private double energy(Phase p) {
        double energy = 0.0;
        for(AtomC a=(AtomC)p.firstAtom(); a!=null; a=a.getNextAtomC()) {
            AtomC nextMoleculeAtom = (AtomC)a.nextMoleculeFirstAtom();
            Potential1 p1 = a.parentMolecule.getP1();
            for(AtomC b=a.getNextAtomC(); b!=nextMoleculeAtom; b=b.getNextAtomC()) {
                energy += p1.getPotential(a,b).energy(a,b);
            }
            Potential2[] p2 = a.parentMolecule.getP2();
            for(AtomC b=nextMoleculeAtom; b!=null; b=b.getNextAtomC()) {
                energy += p2[b.getSpeciesIndex()].getPotential(a,b).energy(a,b);
                if(energy >= Double.MAX_VALUE) {return Double.MAX_VALUE;}
            }
        }
        return energy;
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
        iDisplace = freqDisplace * firstPhase.nMoleculeTotal;
        iVolume = iDisplace + freqVolume;
        iMolecule = iVolume + freqMolecule*firstPhase.nMoleculeTotal;    
        iTotal = iMolecule;
    }
    
    public IntegratorAgent makeAgent(Atom a) {
        return new Agent(a);
    }
    
    private class Agent implements IntegratorAgent {
        public Atom atom;
        public Agent(Atom a) {atom = a;}
    }

}