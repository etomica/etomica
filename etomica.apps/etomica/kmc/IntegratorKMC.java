package etomica.kmc;

import etomica.api.IIntegrator;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.dimer.IntegratorDimerMin;
import etomica.dimer.IntegratorDimerRT;
import etomica.integrator.IntegratorBox;
import etomica.potential.PotentialMaster;
import etomica.space.ISpace;
import etomica.util.Arrays;

public class IntegratorKMC extends IntegratorBox{

    IIntegrator [] integrators;
    IPotentialMaster potentialMaster;
    double temperature;
    private final ISpace space;
    IRandom random;
    ISimulation sim;
    ISpecies [] species;
    
    double[] saddleVib;
    double[] saddleEnergies;
    double[] rates;
    double tau;
    double msd;
    double beta;
    double minEnergy;
    double minVib;
    
    public IntegratorKMC(ISimulation _sim, IPotentialMaster _potentialMaster, double _temperature, IRandom _random, ISpecies [] species, ISpace _space){
        super(_potentialMaster, _temperature);
        this.potentialMaster = _potentialMaster;
        this.temperature = _temperature;
        this.space = _space;
        this.random = _random;
        this.sim = _sim;
        
        // TODO Auto-generated constructor stub
    }

    
    @Override
    protected void doStepInternal() {
        // TODO Auto-generated method stub
        
        /*
        1. Read in energy, vibdata of current state (should be minimum).
        
        2. Offset positions of atoms and spawn DimerSearch simulation X50.
        
        3. Wait....
        
        4. Read in energy, vibdata of dimer searches.  Get activation energy, compute rates, residence time.
        
        5. Choose transition to use, based on rate and probability distribution.
        
        6. Run TWO minimum searches, saving XYZ output, MSD output.
        
        7. Read in energy, vibdata of minimum searches.  
        
        8. Confirm one of the energies is close to current state energy, use opposite minimum, update atom positions to new minimum.
        
        9. Adjust dimer species.
        
        10. Rinse, repeat.
        
         */
        
    }
    
    public void setup(){
        saddleVib = new double[50];
        saddleEnergies = new double[50];
        rates = new double[50];
        
        beta = 1.0/(temperature*1.3806503E-023);
        
        IntegratorDimerMin integratorMinA = new IntegratorDimerMin(sim, potentialMaster, species, true, space);
        IntegratorDimerMin integratorMinB = new IntegratorDimerMin(sim, potentialMaster, species, false, space);
        IntegratorDimerRT integratorDimer = new IntegratorDimerRT(sim, potentialMaster, species, space);
        integrators = (IIntegrator[])Arrays.addObject(integrators,integratorMinA);
        integrators = (IIntegrator[])Arrays.addObject(integrators,integratorMinB);
        integrators = (IIntegrator[])Arrays.addObject(integrators,integratorDimer);
        
        
    }

    public void calcRates(){
        //convert energies to Joules and use hTST
        
        for(int i=0; i<rates.length; i++){
            saddleEnergies[i] = saddleEnergies[i] * 1.60217646E-019;
            rates[i] = (minVib / saddleVib[i]) * Math.exp( -(saddleEnergies[i] - minEnergy)*beta);
        }
    }


}
