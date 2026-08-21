package etomica.virial.simulations.theta;


import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.NeighborIteratorCellFaster;
import etomica.potential.BondingInfo;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManagerCell;
import etomica.species.SpeciesManager;

public class NeighborCellManagerStar extends NeighborCellManager {
    public final AtomType coretype;

    public NeighborCellManagerStar(SpeciesManager sm, Box box, int cellRange, BondingInfo bondingInfo, AtomType coretype) {
        super(sm, box, cellRange, bondingInfo, true);
    this.coretype= coretype;
    }

    @Override
    public NeighborIterator makeNeighborIterator() {
        // return new NeighborIteratorCell(this, bondingInfo, isPureAtoms, box);
        return new NeighborIteratorCellFaster(this, box,coretype);
    }

    @Override
    public void assignCellAll() {
        super.assignCellAll();
        for(IAtom a:box.getLeafList()){
            if(a.getType()==coretype){
                removeAtom(a);
            }
        }
    }

    @Override
    public void updateAtom(IAtom atom) {
        if(atom.getType()!=coretype)
            super.updateAtom(atom);




    }
}
