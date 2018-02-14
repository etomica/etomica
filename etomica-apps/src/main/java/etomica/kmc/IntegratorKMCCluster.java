/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.kmc;

import etomica.action.BoxImposePbc;
import etomica.action.WriteConfiguration;
import etomica.action.XYZWriter;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorBoxDependent;
import etomica.atom.iterator.AtomIteratorLeafFilteredType;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.data.meter.MeterMeanSquareDisplacement;
import etomica.dimer.IntegratorDimerMin;
import etomica.dimer.PotentialMasterListDimer;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorListenerAction;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Joule;
import etomica.util.random.IRandom;

import java.io.*;

public class IntegratorKMCCluster extends IntegratorBox{

    private static final long serialVersionUID = 1L;
    IntegratorDimerMin integratorMin1, integratorMin2;
    PotentialMaster potentialMaster;
    double temperature;
    private final Space space;
    IRandom random;
    Simulation sim;
    ISpecies [] species;
    Vector[] minPosition, currentSaddle, previousSaddle;
    double[] saddleVib;
    double[] saddleEnergies;
    double[] rates;
    double tau;
    Vector msd;
    double beta;
    double massSec;
    double minEnergy;
    double freqProd;
    double saddleEnergy;
    double minVib;
    boolean search;
    int searchNum;
    int kmcStep;
    int totalSearches;
    SimulationGraphic graphic;
    XYZWriter xyzMin1, xyzMin2;
    BoxImposePbc imposePbc;
    MeterMeanSquareDisplacement msd1, msd2;
    FileReader fileReader, writeTau;
    BufferedReader buffReader;
    FileWriter writer;
    
    public IntegratorKMCCluster(Simulation _sim, PotentialMaster _potentialMaster, double _temperature, int _totalSearches, IRandom _random, ISpecies [] _species, Space _space, Box box){
        super(_potentialMaster, _temperature, box);
        
        this.potentialMaster = _potentialMaster;
        this.temperature = _temperature;
        this.space = _space;
        this.random = _random;
        this.sim = _sim;
        this.species = _species;
        this.totalSearches = _totalSearches;

        
                
        // TODO Auto-generated constructor stub
    }

    protected void doStepInternal(){
        loadConfiguration((kmcStep-1)+"");

        for(int i=0; i<totalSearches; i++){
            try {
                FileWriter goWriter = new FileWriter(i+".go");
                goWriter.close();
            } catch (IOException e1) {        }
        }
        
        while(true){
            boolean success = true;
            for(int i=0; i<totalSearches; i++){
                if(!new File(i+".done").exists()){
                    success = false;
                }  
            }
            if(success){
                break;
            }
            try {
                Thread.sleep(20000);
            } catch (InterruptedException e){        }
        }
        for(int i=0; i<totalSearches; i++){
            new File(i+".go").delete();            
        }
        for(int i=0; i<totalSearches; i++){
            new File(i+".done").delete();            
        }
        
        searchNum = 0;
        for(int i=0; i<saddleEnergies.length; i++){

            loadConfiguration("s_"+i+"_saddle");
            try {
                fileReader = new FileReader("s_"+i+"_s_ev");
                buffReader = new BufferedReader(fileReader);
                saddleEnergy = Double.parseDouble(buffReader.readLine());
                freqProd = Double.parseDouble(buffReader.readLine());
            }
            catch(IOException e) {}

            if(checkUniqueSaddle()){
                System.out.println("Good search "+i+", adding saddle data.");
                saddleEnergies[i] = saddleEnergy;
                saddleVib[i] = freqProd;
            }
            searchNum++;
        }
    
        calcRates();
        int rateNum = chooseRate();
        System.out.println("Rate "+rateNum+" is chosen.");
        System.out.println("Step "+(kmcStep-1)+": tau is "+tau);
                
        //Minimum Search with random transition
        integratorMin1.setFileName("s_"+rateNum);
        try {
            integratorMin1.reset();
        } catch (ConfigurationOverlapException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        xyzMin1.setFileName((kmcStep-1)+"_s-A.xyz");
        writeConfiguration((kmcStep-1)+"_s");      
        
        for(int j=0;j<1000;j++){
            integratorMin1.doStep();
            if(integratorMin1.minFound){
                break;
            }
        }
        msdCalc(msd1);
        if(checkMin()){
            minEnergy = integratorMin1.e0;
            minVib = integratorMin1.vib.getProductOfFrequencies();
            writeConfiguration(kmcStep+"");
            writeConfiguration("searchStart");
            setInitialStateConditions(minEnergy, minVib);
            System.out.println("Good minimum found. Computing MSD for other direction...");
            integratorMin2.setFileName("s_"+rateNum);
            try {
                integratorMin2.reset();
            } catch (ConfigurationOverlapException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            xyzMin2.setFileName((kmcStep-1)+"_s-B.xyz");
            for(int j=0;j<1000;j++){
                integratorMin2.doStep();
                if(integratorMin2.minFound){
                    break;
                }
            }
            msdCalc(msd2);
        }else{
            integratorMin2.setFileName("s_"+rateNum);
            //rename minimum 1 search XYZ file
            try {
                integratorMin2.reset();
            } catch (ConfigurationOverlapException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            xyzMin2.setFileName((kmcStep-1)+"_s-A-right.xyz");
            for(int j=0;j<1000;j++){
                integratorMin2.doStep();
                if(integratorMin2.minFound){
                    break;
                }
            }
            msdCalc(msd2);
            minEnergy = integratorMin2.e0;
            minVib = integratorMin2.vib.getProductOfFrequencies();
            writeConfiguration(kmcStep+"");
            writeConfiguration("searchStart");
            setInitialStateConditions(minEnergy, minVib);
            System.out.println("Good minimum found on second attempt.");
        }
                
        try {
            writer = new FileWriter("tau-msd.dat", true);
            writer.write("-step "+kmcStep+"\n"+"tau: "+tau+"\n"+"msd: "+msd+"\n");
            writer.close();
            
            FileWriter writer2 = new FileWriter(kmcStep+"_ev");
            writer2.write(minEnergy+"\n"+minVib);
            writer2.close();
            
        }catch(IOException e) {
            
        }
        clearRatesandEnergies();
        msd1.reset();
        msd2.reset();
        kmcStep++;
    }
    
    public void setup(){
        search = true;
        saddleVib = new double[totalSearches];
        saddleEnergies = new double[totalSearches];
        massSec = Math.sqrt(species[0].getAtomType(0).getMass()) * 0.000000000001;
        msd = space.makeVector();
        tau = 0;
        searchNum = 0;
        kmcStep = 1;
        imposePbc = new BoxImposePbc(box, space);
        rates = new double[totalSearches];
        beta = 1.0/(temperature*1.3806503E-023);
        currentSaddle = new Vector[box.getMoleculeList().getMoleculeCount()];
        previousSaddle = new Vector[box.getMoleculeList().getMoleculeCount()];
        for(int i=0; i<currentSaddle.length; i++){
            currentSaddle[i] = space.makeVector();
            previousSaddle[i] = space.makeVector();
        }
        
        createIntegrators();

    }
    
    public void setInitialStateConditions(double energy, double vibFreq){
        minEnergy = energy;
        minVib = vibFreq;
        
        IMoleculeList loopSet2 = box.getMoleculeList();
        minPosition = new Vector[loopSet2.getMoleculeCount()];
        for(int i=0; i<minPosition.length; i++){
            minPosition[i] = space.makeVector();
        }
        
        for(int i=0; i<loopSet2.getMoleculeCount(); i++){
            minPosition[i].E(loopSet2.getMolecule(i).getChildList().getAtom(0).getPosition());
        }  
    }
    
    public void randomizePositions(){
        Vector workVector = space.makeVector();
        IMoleculeList loopSet3 = box.getMoleculeList(species[0]);
        Vector[] currentPos = new Vector[loopSet3.getMoleculeCount()];
        double offset = 0;
        for(int i=0; i<currentPos.length; i++){
            currentPos[i] = space.makeVector();
            currentPos[i] = (loopSet3.getMolecule(i).getChildList().getAtom(0).getPosition());
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
        minEnergy = Joule.UNIT.fromSim(minEnergy);
        for(int i=0; i<rates.length; i++){
            if(saddleEnergies[i]==0){continue;}
            saddleEnergies[i] = Joule.UNIT.fromSim(saddleEnergies[i]);
            rates[i] = (minVib / saddleVib[i] / massSec) * Math.exp( -(saddleEnergies[i] - minEnergy)*beta);
            rateSum += rates[i];
        }
        //compute residence time
        tau += -Math.log(random.nextDouble())/rateSum;

    }
    
    public void msdCalc(MeterMeanSquareDisplacement msdArray){
        double sum = 0;
        Vector[] msdVect = msdArray.getDataAsArray();
        for(int i=0; i<msdVect.length; i++){
            sum = msd.getX(0);
            sum += msdVect[i].getX(0)*msdVect[i].getX(0);
            msd.setX(0, sum);
            
            sum = msd.getX(1);
            sum += msdVect[i].getX(1)*msdVect[i].getX(1);
            msd.setX(1, sum);
            
            sum = msd.getX(2);
            sum += msdVect[i].getX(2)*msdVect[i].getX(2);
            msd.setX(2, sum);
        }
    }
    
    public void clearRatesandEnergies(){
        for(int i=0; i<rates.length; i++){
            rates[i] = 0.0;
            saddleEnergies[i] = 0.0;
            saddleVib[i] = 0.0;
        }
    }

    public int chooseRate(){
        int rt = 0;
        double sum = 0;
        double rand = random.nextDouble();
        for(int q=0; q<rates.length; q++){
            sum += rates[q];
        }
        double sumgrt = 0;
        for(int i=0; i<rates.length; i++){
            sumgrt += rates[i];
            if(rand*sum<=sumgrt){
                rt = i;
                System.out.println("-----Choosing a rate-----");
                for(int l=0; l<rates.length; l++){ 
                    System.out.println("Rate "+l+": "+rates[l]+", v: "+minVib / saddleVib[i] / massSec);
                }
                System.out.println("Sum:    "+sum);
                System.out.println("-------------------------");
                System.out.println(rand*sum+" <= "+sumgrt);
                break;
            }
        }
        return rt;
    }
    
    private boolean checkUniqueSaddle(){    
        for(int p=0; p<box.getMoleculeList().getMoleculeCount(); p++){
            currentSaddle[p].E(box.getMoleculeList().getMolecule(p).getChildList().getAtom(0).getPosition());
        }
        for(int i=0; i<searchNum; i++){
            double positionDiff = 0;
            loadConfiguration("s_"+i+"_saddle");
            for(int j=0; j<box.getMoleculeList().getMoleculeCount(); j++){
                previousSaddle[j].E(box.getMoleculeList().getMolecule(j).getChildList().getAtom(0).getPosition());
                previousSaddle[j].ME(currentSaddle[j]);
                positionDiff += previousSaddle[j].squared();
            }
            if(positionDiff < 0.5){
                System.out.println("Duplicate saddle found.");
                return false;
            }
        }
        System.out.println("Unique saddle found.");
        return true;  
    }
    
    public boolean checkMin(){
        Vector workVector = space.makeVector();
        double positionDiff=0;
        for(int i=0; i<box.getMoleculeList().getMoleculeCount(); i++){
            workVector.Ev1Mv2(minPosition[i],box.getMoleculeList().getMolecule(i).getChildList().getAtom(0).getPosition());
            positionDiff += workVector.squared();
        }
        if(positionDiff > 0.5){return true;}
        return false;
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
    
    public void createIntegrators() {
        integratorMin1 = new IntegratorDimerMin(sim, potentialMaster, species, true, space, box);
        integratorMin2 = new IntegratorDimerMin(sim, potentialMaster, species, false, space, box);

        if (potentialMaster instanceof PotentialMasterListDimer) {
            integratorMin2.getEventManager().addListener(((PotentialMasterList) potentialMaster).getNeighborManager(box));
        }

        xyzMin1 = new XYZWriter(box);
        xyzMin2 = new XYZWriter(box);
        xyzMin1.setIsAppend(true);
        xyzMin2.setIsAppend(true);

        IntegratorListenerAction xyzMin1Listener = new IntegratorListenerAction(xyzMin1);
        xyzMin1Listener.setInterval(5);
        IntegratorListenerAction xyzMin2Listener = new IntegratorListenerAction(xyzMin2);
        xyzMin2Listener.setInterval(5);
        integratorMin1.getEventManager().addListener(xyzMin1Listener);
        integratorMin2.getEventManager().addListener(xyzMin2Listener);
        IntegratorListenerAction imposePbc1Listener = new IntegratorListenerAction(imposePbc);
        imposePbc1Listener.setInterval(1);
        IntegratorListenerAction imposePbc2Listener = new IntegratorListenerAction(imposePbc);
        imposePbc2Listener.setInterval(1);
        integratorMin1.getEventManager().addListener(imposePbc1Listener);
        integratorMin2.getEventManager().addListener(imposePbc2Listener);

        //Limit MSD calculation to a specific species
        AtomIterator aif = new AtomIteratorLeafFilteredType(box, species[0].getAtomType(0));
        msd1 = new MeterMeanSquareDisplacement(space, integratorMin1);
        msd2 = new MeterMeanSquareDisplacement(space, integratorMin2);
        msd1.setIterator((AtomIteratorBoxDependent) aif);
        msd2.setIterator((AtomIteratorBoxDependent) aif);
    }
    
    
}
