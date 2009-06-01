package etomica.api;

public interface IBoxMoleculeCountEvent extends IBoxEvent {

    public ISpecies getSpecies();
    
    public int getCount();
}
