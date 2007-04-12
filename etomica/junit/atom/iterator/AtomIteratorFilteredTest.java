package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomFilter;
import etomica.atom.AtomLeaf;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorFiltered;
import etomica.space3d.Space3D;

/**
 * Unit test of AtomIteratorFiltered
 * 
 * @author David Kofke
 *  
 */
public class AtomIteratorFilteredTest extends IteratorTestAbstract {

    /**
     * Filter that rejects every n-th atom, knowing that ordinals are numbered 1
     * to n
     */
    public class MyFilter implements AtomFilter {

        int n = 0;

        public boolean accept(IAtom a) {
            if(n == 0) return true;
            return ((a.getIndex()+1) % n != 0);
        }

    }

    public void testIterator() {

        AtomArrayList list = new AtomArrayList();

        AtomIteratorArrayListSimple listIterator = new AtomIteratorArrayListSimple(list);

        MyFilter filter = new MyFilter();

        AtomIteratorFiltered iterator = AtomIteratorFiltered.makeIterator(
                listIterator, filter);

        for (filter.n = 0; filter.n < 10; filter.n++) {
            list.clear();
            for (int j = 0; j < 20; j++) {
                LinkedList iterateList = generalIteratorMethodTests(iterator);
                int count = (filter.n == 0) ? j : (filter.n - 1)*j; 
                //System.out.println(list.size()+ " " + iterateList.size() + " " + count);
                assertEquals(iterateList.size(), count);
                //add n atoms at a time, numbering ordinals 1 to n
                int n = Math.max(filter.n,1);
                for (int k = 0; k < n; k++) {
                    IAtom atom = new AtomLeaf(Space3D.getInstance());
                    atom.setIndex(k);
                    list.add(atom);
                }
            }
        }
    }

}
