package etomica;

public abstract class Potential1Soft extends Potential1 {
    
    public abstract Space.Vector gradient(Atom atom);
    
    public Potential1Soft(Simulation sim) {super(sim);}
    
    public void calculate(IteratorDirective id, PotentialCalculation pc) {
        if( !(pc instanceof Potential1Calculation) ) return;
        iterator.reset(id);
        ((Potential1Calculation)pc).calculate(iterator, this); 
    }
}