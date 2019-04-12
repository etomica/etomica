package etomica.virial.simulations.KnottedPolymer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialMolecular;
import etomica.space.Space;

import java.util.HashMap;


public class P2HSPolymer extends PotentialMolecular {
    final IPotentialAtomic potential;
    private final HashMap<Box, Api1ACell> iteratorMap;
    private final HashMap<Box, NeighborCellManager> boxMapper;
    public double range = 1.00; // consider the atom diameter to be 1.00
    private final double r2c = range * range;
    private boolean atomOutBox = false;
    //private double[][] energyMoleculeArray;

    public P2HSPolymer(int nBody, Space space, IPotentialAtomic potential) {
        super(nBody, space);
        this.potential = potential;
        this.boxMapper = new HashMap<>();
        this.iteratorMap = new HashMap<>();
    }

    @Override
    public double getRange() {
        return range;
    }

    @Override
    public double energy(IMoleculeList molecules) {
        if (!boxMapper.keySet().isEmpty() && !iteratorMap.keySet().isEmpty()) {
            Box box = boxMapper.keySet().iterator().next();
            Api1ACell api1ACell = iteratorMap.get(box);

//            energyMoleculeArray = new double[box.getMoleculeList().size()][box.getMoleculeList().size()];
//
//            for (int i=0; i<box.getMoleculeList().size(); i++){
//                for (int j=0; j<box.getMoleculeList().size(); j++){
//                    energyMoleculeArray[i][j] = 0.0;
//                }
//            }

            for (IAtomList ial = api1ACell.next(); ial != null; ial = api1ACell.next()) {
                IAtom atom1 = ial.get(0);
                IAtom atom2 = ial.get(1);

                if (atom1.getParentGroup() == atom2.getParentGroup()) {
                    continue;
                }
//                if (energyMoleculeArray[molecule1.getIndex()][molecule2.getIndex()] == Double.POSITIVE_INFINITY){
////                    return Double.POSITIVE_INFINITY;
////                }

                double r2 = atom1.getPosition().Mv1Squared(atom2.getPosition());
                if (r2 < r2c) {
//                    energyMoleculeArray[molecule2.getIndex()][molecule1.getIndex()] = Double.POSITIVE_INFINITY;
//                    energyMoleculeArray[molecule1.getIndex()][molecule2.getIndex()] = Double.POSITIVE_INFINITY;
                    return Double.POSITIVE_INFINITY;
                }

            }

        }
//        return energyMoleculeArray[molecules.get(0).getIndex()][molecules.get(1).getIndex()];
        return 0.0;
    }

    @Override
    public void setBox(Box box) {
        if (!box.getLeafList().getAtoms().isEmpty()) {
            for (int i = 0; i < box.getLeafList().size(); i++) {
                IAtom atom = box.getLeafList().get(i);
                if (!box.getBoundary().getShape().contains(atom.getPosition())) {
                    atomOutBox = true;
                    break;
                }
            }
        }

        if (!atomOutBox) {
            Api1ACell api1ACell = iteratorMap.get(box);
            NeighborCellManager neighborCellManager = boxMapper.get(box);
            if (neighborCellManager == null) {
                neighborCellManager = new NeighborCellManager(box, 8.0);
                boxMapper.put(box, neighborCellManager);
                api1ACell = new Api1ACell(getRange(), box, neighborCellManager);
                iteratorMap.put(box, api1ACell);
            }
//            neighborCellManager.setCellRange(1);
            neighborCellManager.assignCellAll();
            api1ACell.setTarget(box.getMoleculeList().get(0).getChildList().get(0));
            api1ACell.reset();
        }
    }
}
