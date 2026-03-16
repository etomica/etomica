package etomica.nbr.cell;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.BondingInfo;
import etomica.potential.compute.NeighborIterator;
import etomica.space.Space;
import etomica.space.Vector;

import java.util.Arrays;

public class NeighborIteratorCellMulti implements NeighborIterator {

    private final NeighborCellManagerMulti cellManager;
    private final BondingInfo bondingInfo;

    private final Box box;
    private final Space space;



    public NeighborIteratorCellMulti(NeighborCellManagerMulti cellManager, BondingInfo bondingInfo,  Box box) {
        this.cellManager = cellManager;
        this.bondingInfo = bondingInfo;

        this.box = box;
        this.space = box.getSpace();

    }
    public static boolean doDebug ;



    @Override
    public void iterUpNeighbors(int iAtom, NeighborConsumer consumer) {

        IAtomList atoms = box.getLeafList();

        int[] atomCell = cellManager.getAtomCell();
        int[] cellNextAtom = cellManager.getCellNextAtom();
        int[] cellOffsets = cellManager.getCellOffsets();
        int[] cellCoordinate = cellManager.atomCellCoordinate;

        int[][] cellLastAtom = cellManager.getCellLastAtom();
        int cellRange = cellManager.cellRange;
        int [] numCells = cellManager.numCells;
        int totalCells= cellLastAtom[0].length;
        IAtom atom1 = atoms.get(iAtom);
        Vector rij = space.makeVector();
        int m = atom1.getParentGroup().getIndex();
        int [][] offset = cellManager.moleculeCellOffsets;
        int debugAtom1= 130;
        int debugAtom2= 1119;
        int yoffset=128- ((int) ((0.5 ) * (numCells[0] - 2 * cellRange)));


        if(doDebug && iAtom==debugAtom1){
            int c = cellCoordinate[debugAtom2];
            int []cxyz= new int[3];

            cxyz[0] = c%256-yoffset;
            cxyz[1]= (c/256)%256-yoffset;
            cxyz[2]=( (c/256)/256)%256-yoffset;
            IAtom atomdebugAtom1 = atoms.get(debugAtom1);
            IAtom atomdebugAtom2 = atoms.get(debugAtom2);
            IAtom atom0 = atoms.get(0);
            IAtom atom641= atoms.get(641);


            System.out.println(atomCell[debugAtom2]+" "+cellCoordinate[debugAtom2]+" "+cxyz[0]+" "+cxyz[1]+" "+cxyz[2]);
             c = cellCoordinate[debugAtom1];

            cxyz[0] = c%256-yoffset;
            cxyz[1]= (c/256)%256-yoffset;
            cxyz[2]=( (c/256)/256)%256-yoffset;
            System.out.println(atomCell[debugAtom1]+" "+cellCoordinate[debugAtom1]+" "+cxyz[0]+" "+cxyz[1]+" "+cxyz[2]);
            System.out.println("Atom debugAtom2:   pos=" + atomdebugAtom2.getPosition());
            System.out.println("Atom debugAtom1:   pos=" + atomdebugAtom1.getPosition() );
            System.out.println("Atom 0:   pos=" + atom0.getPosition() );
            System.out.println("Atom 641:   pos=" + atom641.getPosition() );


            System.out.println("CellOffsets: offsets" + Arrays.toString(offset[0]));
            System.out.println("CellOffsets: offsets" + Arrays.toString(offset[1]));
        }

        for (int k= m;k<box.getMoleculeList().size();k++) {
            int[] kOffset= new int[3];
            for(int l =0;l<3;l++){
                kOffset[l]= offset[k][l]-offset[m][l];
            }
            int iCell = atomCell[iAtom];
            int [] cxyz= new int[3];
            int [] minBigCell = new int[]{-1,-1,-1};
            int [] maxBigCell = new int[]{1,1,1};
            if(k!=m){
                int c = cellCoordinate[iAtom];

               cxyz[0] = c%256-kOffset[0]-yoffset;
               cxyz[1]= (c/256)%256-kOffset[1]-yoffset;
               cxyz[2]=( (c/256)/256)%256-kOffset[2]-yoffset;
                if (doDebug && iAtom==debugAtom1){
                    System.out.println("c "+c+" "+cxyz[0]+" "+cxyz[1]+" "+cxyz[2]);
                }
               int bigCellNum = 0;
               int bigCellJump = 1;
               iCell = 0;
               int [] jump= cellManager.jump;
               for(int l =0;l<3;l++) {
                   int y = cxyz[l];

                   if (y < cellRange || y >= numCells[l] - cellRange) {

                       bigCellNum += (y < cellRange ? -1 : 1) * bigCellJump;
                   } else {
                       //if (y == numCells[l] - cellRange) y--;
                       //else if (y == cellRange - 1) y++;
                       iCell += y * jump[l];
                   }
                   bigCellJump *= 5;
               }
                if (bigCellNum!=0) {
                    for(int l =0;l<3;l++) {
                        int y = cxyz[l];
                        if(y>numCells[l]){
                            minBigCell[l]=1;
                        } else if (y>=2*cellRange) {
                            minBigCell[l]=0;

                        }
                        if(y<0){
                            maxBigCell[l]=-1;

                        } else if (y<numCells[l]-2*cellRange) {
                            maxBigCell[l]=0;

                        }
                    }
                    iCell = cellLastAtom[k].length - 63 + bigCellNum;
                }
            }
            if (doDebug && iAtom==debugAtom1){
                System.out.println(k+" iCell "+iCell);
            }

            int jStart= k ==m? cellNextAtom[iAtom]:cellLastAtom[k][iCell];

            for (int j = jStart; j > -1; j = cellNextAtom[j]) {
                IAtom atom2 = atoms.get(j);
                boolean skipIntra = bondingInfo.skipBondedPair(false, atom1, atom2);
                if (skipIntra) continue;
                int n = bondingInfo.n(false, atom1, atom2);

                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());


                consumer.accept(atom2, rij, n);
            }


            if (iCell >= totalCells - 125) {

                /*for (int x = -1; x <= 1; x++) {
                    for (int y = -1; y <= 1; y++) {
                        for (int z = -1; z <= 1; z++) {

                            int jCell = iCell + z * 25 + y * 5 + x;*/
                if(doDebug && iAtom==debugAtom1){
                    System.out.println(minBigCell[0]+" "+maxBigCell[0]+" "+minBigCell[1]+" "+maxBigCell[1]+" "+minBigCell[2]+" "+maxBigCell[2]);
                }
                for(int x = minBigCell[0];x<= maxBigCell[0];x++){
                    for(int y = minBigCell[1];y<= maxBigCell[1];y++){
                        for(int z = minBigCell[2];z<= maxBigCell[2];z++) {
                            int jCell=totalCells-63+z*25+y*5+x;
                            if(doDebug && iAtom==debugAtom1) {
                                System.out.println("j " + jCell);
                            }

                            if (jCell==iCell || (jCell < iCell && k==m )) continue;

                            for (int j = cellLastAtom[k][jCell]; j > -1; j = cellNextAtom[j]) {
                                IAtom atom2 = atoms.get(j);


                                //Vector rij = space.makeVector();
                                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());


                                consumer.accept(atom2, rij, 0);

                            }
                        }
                    }
                }
                int c = iCell - totalCells + 125;
                int x = c % 5;
                int y = ((c - x) / 5) % 5;
                int z = (((c - x) / 5) - y) / 5;
                int cr = cellManager.cellRange;
                int minX = cr, minY = cr, minZ = cr;
                int maxX = cellManager.numCells[0] - cr-1, maxY = cellManager.numCells[1] - cr-1, maxZ = cellManager.numCells[2] - cr-1;
                if(k==m) {
                    if (x == 1) maxX = minX + cr - 1;
                    else if (x == 3) minX = maxX - cr + 1;
                    if (y == 1) maxY = minY + cr - 1;
                    else if (y == 3) minY = maxY - cr + 1;
                    if (z == 1) maxZ = minZ + cr - 1;
                    else if (z == 3) minZ = maxZ - cr + 1;
                }
                else {
                    minX=Math.max(minX,cxyz[0]-2);
                    maxX=Math.min(maxX,cxyz[0]+2);
                    minY=Math.max(minY,cxyz[1]-2);
                    maxY=Math.min(maxY,cxyz[1]+2);
                    minZ=Math.max(minZ,cxyz[2]-2);
                    maxZ=Math.min(maxZ,cxyz[2]+2);
                }
                for (int cX = minX; cX <= maxX; cX++) {
                    for (int cY = minY; cY <= maxY; cY++) {
                        for (int cZ = minZ; cZ <= maxZ; cZ++) {
                            int jCell = cX + cY * cellManager.numCells[0] + cZ * cellManager.numCells[1] * cellManager.numCells[0];
                            if(doDebug && iAtom==debugAtom1 && k!=m){
                                System.out.println(cX+" "+cY+" "+cZ+" "+jCell);
                            }
                            for (int j = cellLastAtom[k][jCell]; j > -1; j = cellNextAtom[j]) {
                                IAtom atom2 = atoms.get(j);


                                // Vector rij = space.makeVector();
                                rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());


                                consumer.accept(atom2, rij, 0);

                            }
                        }
                    }
                }
            } else {
                for (int ico = 0; ico < cellManager.numCellOffsets; ico++) {
                    int cellOffset = cellOffsets[ico];
                    int jCell = iCell + cellOffset;
                    if(jCell>= cellLastAtom[k].length){
                        System.out.println(iCell+ " "+ ico+" "+ jCell+" "+cellLastAtom[k].length );

                    }
                    for (int j = cellLastAtom[k][jCell]; j > -1; j = cellNextAtom[j]) {
                        IAtom atom2 = atoms.get(j);
                        rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                        consumer.accept(atom2, rij, 0);
                    }
                    if(k!=m) {
                        jCell = iCell - cellOffset;
                        for (int j = cellLastAtom[k][jCell]; j > -1; j = cellNextAtom[j]) {
                            IAtom atom2 = atoms.get(j);
                            rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                            consumer.accept(atom2, rij, 0);
                        }
                    }
                }
                if(k!=m) {
                    int minX = -1, minY = -1, minZ = -1;
                    int maxX = 1, maxY = 1, maxZ = 1;
                    if(cxyz[0]>cellRange+1) minX=0;
                    if(cxyz[0]<numCells[0]-2*cellRange) maxX=0;
                    if(cxyz[1]>cellRange+1) minY=0;
                    if(cxyz[1]<numCells[1]-2*cellRange) maxY=0;
                    if(cxyz[2]>cellRange+1) minZ=0;
                    if(cxyz[2]<numCells[2]-2*cellRange) maxZ=0;

                    for( int iX=minX;iX<=maxX;iX++){
                        for( int iY=minY;iY<=maxY;iY++){
                            for( int iZ=minZ;iZ<=maxZ;iZ++) {
                                if(iX==0 && iY==0 && iZ==0) continue;

                                int jCell = totalCells - 63+ iZ*25 +iY*5 +iX;
                                for (int j = cellLastAtom[k][jCell]; j > -1; j = cellNextAtom[j]) {
                                    IAtom atom2 = atoms.get(j);
                                    rij.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                                    consumer.accept(atom2, rij, 0);
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    @Override
    public void iterDownNeighbors(int iAtom, NeighborConsumer consumer) {

    }

    @Override
    public double iterAndSumAllNeighbors(IAtom atom1, SuperNbrConsumer consumer) {
        return 0;
    }

    @Override
    public void iterAllNeighbors(int iAtom, NeighborConsumer consumer) {

    }

}
