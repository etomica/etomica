package simulate;

public interface LatticeSite {
    public AtomP atom();
    public void putAtom(AtomP a);
    public LatticeSite[] neighbors();
}