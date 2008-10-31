package etomica.kmc;

import etomica.action.WriteConfiguration;
import etomica.action.XYZWriter;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.config.ConfigurationFile;
import etomica.dimer.IntegratorDimerMin;
import etomica.dimer.IntegratorDimerRT;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorBox;
import etomica.space.ISpace;

public class IntegratorKMC extends IntegratorBox{

    private static final long serialVersionUID = 1L;
    IntegratorDimerRT integratorDimer;
    IntegratorDimerMin integratorMin1, integratorMin2;
    IPotentialMaster potentialMaster;
    double temperature;
    private final ISpace space;
    IRandom random;
    ISimulation sim;
    ISpecies [] species;
    IVector [] minPosition;
    double[] saddleVib;
    double[] saddleEnergies;
    double[] rates;
    double tau;
    double msd;
    double beta;
    double minEnergy;
    double minVib;
    boolean search;
    int goodSearch;
    int stepCounter;
    int searchlimit;
    SimulationGraphic graphic;
    XYZWriter xyzfile;
    
    public IntegratorKMC(ISimulation _sim, IPotentialMaster _potentialMaster, double _temperature, IRandom _random, ISpecies [] _species, ISpace _space){
        super(_potentialMaster, _temperature);
        this.potentialMaster = _potentialMaster;
        this.temperature = _temperature;
        this.space = _space;
        this.random = _random;
        this.sim = _sim;
        this.species = _species;
        
        searchlimit = 5;
        tau = 0;
        
                
        // TODO Auto-generated constructor stub
    }
 
    @Override
    protected void doStepInternal(){
        
        //Dimer Searches from minimum
        goodSearch = 0;
        while(search){
            loadConfiguration(stepCounter+"");
            randomizePositions();
            try {
                System.out.println("Initializing dimer.");
                integratorDimer.initialize();
            } catch (ConfigurationOverlapException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            integratorDimer.setFileName("s_"+goodSearch);
            System.out.println("Searching...");
            for(int j=0;j<800;j++){
                integratorDimer.doStep();
                if(integratorDimer.saddleFound){
                    System.out.println("Good search "+goodSearch+", adding saddle data.");
                    saddleEnergies[goodSearch] = integratorDimer.saddleEnergy;
                    saddleVib[goodSearch] = integratorDimer.vib.getProductOfFrequencies();
                    goodSearch++;
                    break;
                }
                
            }
            if(goodSearch>searchlimit-1){search = false;}
        }
        search = true;
                
        calcRates();
        int rateNum = chooseRate();
        System.out.println("Rate "+rateNum+" is chosen.");
        
        //Minimum Search with random transition
        integratorMin1.setFileName("s_"+rateNum);
        try {
            integratorMin1.initialize();
        } catch (ConfigurationOverlapException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        integratorMin1.initializeDimer();
        writeConfiguration(stepCounter+"_saddle");
        xyzfile.actionPerformed();
        
        stepCounter++;
        
        for(int j=0;j<1000;j++){
            integratorMin1.doStep();
            if(integratorMin1.minFound){
                break;
            }
        }
        boolean check = checkMin();
        if(check==true){
            minEnergy = integratorMin1.e0;
            minVib = integratorMin1.vib.getProductOfFrequencies();
            writeConfiguration(stepCounter+"");
            setInitialStateConditions(minEnergy, minVib);
            System.out.println("Good minimum found.");
        }else{
            integratorMin2.setFileName("s_"+rateNum);
            try {
                integratorMin2.initialize();
            } catch (ConfigurationOverlapException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            integratorMin2.initializeDimer();
            for(int j=0;j<1000;j++){
                integratorMin2.doStep();
                if(integratorMin2.minFound){
                    break;
                }
            }
            minEnergy = integratorMin1.e0;
            minVib = integratorMin1.vib.getProductOfFrequencies();
            writeConfiguration(stepCounter+"");
            setInitialStateConditions(minEnergy, minVib);
            System.out.println("Good minimum found on second attempt.");
        }
        xyzfile.actionPerformed();
        
    }
    
    public void setup(){
        search = true;
        saddleVib = new double[searchlimit];
        saddleEnergies = new double[searchlimit];
        
        rates = new double[searchlimit];
        beta = 1.0/(temperature*1.3806503E-023);
        stepCounter = 0;              
    }

    public void setInitialStateConditions(double energy, double vibFreq){
        minEnergy = energy;
        minVib = vibFreq;
        
        IAtomSet loopSet2 = box.getMoleculeList();
        minPosition = new IVector[loopSet2.getAtomCount()];
        for(int i=0; i<minPosition.length; i++){
            minPosition[i] = space.makeVector();
        }
        
        for(int i=0; i<loopSet2.getAtomCount(); i++){
            minPosition[i].E(((IAtomPositioned)((IMolecule)loopSet2.getAtom(i)).getChildList().getAtom(0)).getPosition());
        }  
    }
    
    public void setSearchLimit(int limit){
        searchlimit = limit;
    }
    
    public void randomizePositions(){
        IVector workVector = space.makeVector();
        IAtomSet loopSet3 = box.getMoleculeList(species[0]);
        IVector [] currentPos = new IVector [loopSet3.getAtomCount()];
        double offset = 0;
        for(int i=0; i<currentPos.length; i++){
            currentPos[i] = space.makeVector();
            currentPos[i] = (((IAtomPositioned)((IMolecule)loopSet3.getAtom(i)).getChildList().getAtom(0)).getPosition());
            for(int j=0; j<3; j++){
                offset = random.nextGaussian()/10.0;
                if(Math.abs(offset)>0.1){offset=0.1;}
                workVector.setX(j,offset);
            }
            currentPos[i].PE(workVector);
        }
    }
        
    public void calcRates(){
        //convert energies to Joules and use hTST
        double rateSum = 0;
        double rate = 0;
        minEnergy = minEnergy * 1.60217646E-019;
        for(int i=0; i<rates.length; i++){
            saddleEnergies[i] = saddleEnergies[i] * 1.60217646E-019;
            rates[i] = (minVib / saddleVib[i]) * Math.exp( -(saddleEnergies[i] - minEnergy)*beta);
            rateSum += rates[i];
        }
        
        //compute residence time
        tau += -Math.log(random.nextDouble())/rateSum;

    }
    
    public int chooseRate(){
        int rt = 0;
        double sum = 0;
        double rand = random.nextDouble();
        for(int q=0; q<rates.length; q++){
            sum += rates[q];
        }
        for(int i=0; i<rates.length; i++){
            double sumless = 0;
            double sumgreater = 0;
            for(int j=0; j<i+1; j++){
                sumgreater += rates[j];
            }
            for(int k=0; k<i; k++){
                sumless += rates[k];
            }
            if(rand*sum>sumless && rand*sum<=sumgreater){
                rt = i;
                System.out.println("-----Choosing a rate-----");
                for(int l=0; l<rates.length; l++){ 
                    System.out.println("Rate "+l+": "+rates[l]);
                }
                System.out.println("Sum:    "+sum);
                System.out.println("-------------------------");
                System.out.println(sumless+" < "+rand*sum+" <= "+sumgreater);
                
                break;
            }
        }
        return rt;
    }

    public boolean checkMin(){
        boolean goodMin = false;
        IVector workVector = space.makeVector();
        double positionDiff=0;
        for(int i=0; i<box.getMoleculeList().getAtomCount(); i++){
            workVector.Ev1Mv2(minPosition[i],((IAtomPositioned)((IMolecule)box.getMoleculeList().getAtom(i)).getChildList().getAtom(0)).getPosition());
            positionDiff += workVector.squared();
        }
        if(positionDiff > 0.5){goodMin = true;}
        return goodMin;
    }
        
    public void writeConfiguration(String file){
        WriteConfiguration writer = new WriteConfiguration(space);
        writer.setBox(box);
        writer.setConfName(file);
        writer.actionPerformed();
    }
    
    public void loadConfiguration(String file){
        ConfigurationFile config = new ConfigurationFile(file);
        config.initializeCoordinates(box);
    }

    public void createIntegrators(){
        integratorMin1 = new IntegratorDimerMin(sim, potentialMaster, species, true, space);
        integratorMin2= new IntegratorDimerMin(sim, potentialMaster, species, false, space);
        integratorDimer = new IntegratorDimerRT(sim, potentialMaster, species, space);
        
        integratorMin1.setBox(box);
        integratorMin2.setBox(box);
        integratorDimer.setBox(box);
        integratorDimer.setRotNum(0);
        integratorDimer.setOrtho(false, false);
        
        try {
            integratorDimer.initialize();
            integratorMin1.initialize();
            integratorMin2.initialize();
        } catch (ConfigurationOverlapException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        xyzfile = new XYZWriter(box);
        xyzfile.setIsAppend(true);
        xyzfile.setFileName("kmc-lj-2");
    }
    

}