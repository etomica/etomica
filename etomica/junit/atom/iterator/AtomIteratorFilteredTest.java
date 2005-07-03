package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.Atom;
import etomica.Space;
import etomica.atom.AtomFilter;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorFiltered;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.space1d.Space1D;

/**
 * Unit test of AtomIteratorFiltered
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jul 2, 2005 by kofke
 */
public class AtomIteratorFilteredTest extends IteratorTest {

    /**
     * Filter that rejects every n-th atom, knowing that ordinals are numbered 1
     * to n
     */
    public class MyFilter implements AtomFilter {

        int n = 0;

        public boolean accept(Atom a) {
            if(n == 0) return true;
            return (a.node.getOrdinal() % n != 0);
        }

    }

    public void testIterator() {

        Space space = new Space1D();

        AtomList list = new AtomList();

        AtomIteratorListSimple listIterator = new AtomIteratorListSimple(list);

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
                    Atom atom = new Atom(space);
                    atom.node.setOrdinal(k + 1);
                    list.add(atom);
                }
            }
        }
    }

}
