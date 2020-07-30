package etomica.action.activity;

import etomica.action.IAction;
import etomica.action.activity.Controller2.ActivityHandle;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.*;

class ControllerTest {

    @Test
    public void testActivityBasic() {
        Controller2 controller = new Controller2();

        TestActivity act = new TestActivity();
        ActivityHandle handle = controller.addActivity(act, 20, 1);

        assertFalse(act.didPre);
        assertFalse(act.didPost);
        assertEquals(0, act.runCount);

        controller.start();

        try {
            handle.future.get();
        } catch (InterruptedException | ExecutionException e) {
            fail(e);
        }

        assertTrue(act.didPre);
        assertTrue(act.didPost);
        assertEquals(20, act.runCount);
    }

    @Test
    void testActions() {
        Controller2 controller = new Controller2();

        TestActivity act = new TestActivity();
        act.time = 1;
        ActivityHandle handle = controller.addActivity(act, 100, 1);
        List<TestAction> actions = Stream.generate(TestAction::new).limit(50).collect(Collectors.toList());
        controller.start();
        List<CompletableFuture<Void>> actionFutures = actions.stream().map(controller::submitActionInterrupt).collect(Collectors.toList());
        try {
            handle.future.get();
        } catch (InterruptedException | ExecutionException e) {
            fail(e);
        }
        actions.forEach(a -> assertTrue(a.ran));
        actionFutures.forEach(f -> f.isDone());
    }

    static class TestActivity implements Activity2 {
        boolean didPre = false;
        boolean didPost = false;
        int runCount = 0;
        int time = 0;

        @Override
        public void preAction() {
            didPre = true;
        }

        @Override
        public void postAction() {
            didPost = true;
        }

        @Override
        public void actionPerformed() {
            try {
                Thread.sleep(time);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            runCount++;
        }
    }

    static class TestAction implements IAction {
        boolean ran = false;
        @Override
        public void actionPerformed() {
            ran = true;
        }
    }
}