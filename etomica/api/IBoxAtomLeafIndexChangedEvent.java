package etomica.api;

public interface IBoxAtomLeafIndexChangedEvent extends IBoxAtomEvent {

    public abstract int getOldIndex();

}