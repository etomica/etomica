package etomica.action.activity;

import etomica.action.IAction;
import etomica.action.controller.Controller;
import etomica.action.controller.Controller.ActivityHandle;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.*;

@Timeout(value = 1, unit = TimeUnit.SECONDS)
class ControllerTest {

    private static void sleep(long millis) {
        try {
            Thread.sleep(millis);
        } catch (InterruptedException e) {
            fail(e);
        }
    }

    @Test
    public void testActivityBasic() {
        Controller controller = new Controller();

        TestActivity act = new TestActivity();
        ActivityHandle handle = controller.addActivity(act, 20, 1);

        assertFalse(act.didPre);
        assertFalse(act.didPost);
        assertEquals(0, act.runCount);

        controller.start();
        controller.unpause();

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
        Controller controller = new Controller();

        TestActivity act = new TestActivity();
        act.time = 1;
        ActivityHandle handle = controller.addActivity(act, 100, 1);
        List<TestAction> actions = Stream.generate(TestAction::new).limit(50).collect(Collectors.toList());
        controller.start();
        controller.unpause();
        List<CompletableFuture<Void>> actionFutures = actions.stream().map(controller::submitActionInterrupt).collect(Collectors.toList());
        try {
            handle.future.get();
        } catch (InterruptedException | ExecutionException e) {
            fail(e);
        }
        actions.forEach(a -> assertTrue(a.ran));
        actionFutures.forEach(f -> f.isDone());
    }

    @Test
    @DisplayName("Actions submitted from within an action should execute before the next activity step")
    void testSubmitActionInAction() {
        Controller controller = new Controller();

        TestActivity act = new TestActivity();
        act.time = 1;

        ActivityHandle handle = controller.addActivity(act, 100, 0);

        controller.start();
        controller.unpause();

        sleep(2);

        controller.submitActionInterrupt(() -> {
            int curStep = act.runCount;

            for (int i = 0; i < 100; i++) {
                controller.submitActionInterrupt(() -> {
                    assertEquals(curStep, act.runCount);
                });
            }
        });

        AtomicInteger step = new AtomicInteger(-1);
        controller.submitActionInterrupt(() -> {
            step.set(act.runCount);
        }).whenComplete((res, ex) -> {
            for (int i = 0; i < 100; i++) {
                controller.submitActionInterrupt(() -> {
                    assertEquals(step.get(), act.runCount);
                });
            }
        });

    }

    static class TestActivity implements Activity {
        volatile boolean didPre = false;
        volatile boolean didPost = false;
        volatile int runCount = 0;
        volatile int time = 0;

        @Override
        public void preAction() {
            didPre = true;
        }

        @Override
        public void postAction() {
            didPost = true;
        }

        @Override
        public void restart() {

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