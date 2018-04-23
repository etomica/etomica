package etomica.action;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Created by kofke on 11/11/17.
 */
class ActionGroupSeriesTest {

    private ActionGroupSeries ags = new ActionGroupSeries();
    private int sum = 0;
    private IAction a1, a2, zero;

    @BeforeEach
    void setUp() {
        a1 = new Add1();
        a2 = new Add2();
        zero = new Zero();
    }

    @Test
    void testActionPerformed() {
        ags.actionPerformed();
        assertEquals(0, sum);

        ags.addAction(a1);
        ags.actionPerformed();
        assertEquals(1,sum);

        sum = 0;
        ags.addAction(a2);
        ags.actionPerformed();
        assertEquals(3,sum);

        sum = 0;
        ags.addAction(a1);
        ags.actionPerformed();
        assertEquals(4,sum);

        sum = 0;
        ags.removeAction(a1);
        ags.actionPerformed();
        assertEquals(3,sum);

        sum = 0;
        ags.removeAction(a1);
        ags.actionPerformed();
        assertEquals(2,sum);

        sum = 0;
        assertTrue(ags.removeAction(a2));
        ags.actionPerformed();
        assertEquals(0,sum);

        //test return false for removing action not in list
        assertFalse(ags.removeAction(a2));

        //test constructor
        sum = 0;
        ags = new ActionGroupSeries(new IAction[] {a2, a2, a2, a2, a1, a2, a2});
        ags.actionPerformed();
        assertEquals(13,sum);

        //test getAllActions
        IAction[] actions = ags.getAllActions();
        assertEquals(7, actions.length);

        //test that remove deletes first instance, not last
        sum = 0;
        ags = new ActionGroupSeries(new IAction[] {a2,zero,a2});
        ags.removeAction(a2);
        ags.actionPerformed();
        assertEquals(2, sum);

    }

    private class Add1 implements IAction {
        public void actionPerformed() {sum++;}
    }
    private class Add2 implements IAction {
        public void actionPerformed() {sum += 2;}
    }
    private class Zero implements IAction {
        public void actionPerformed() {sum = 0;}
    }

}
