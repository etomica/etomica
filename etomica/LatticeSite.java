package simulate;

public interface LatticeSite {
    public AtomL atom();
    public void putAtom(AtomL a);
    public LatticeSite[] neighbors();
}