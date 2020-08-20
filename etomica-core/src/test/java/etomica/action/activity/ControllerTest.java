package etomica.action.activity;

import etomica.action.IAction;
import etomica.action.controller.Activity;
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
        act.steps = 20;
        controller.setSleepPeriod(1);
        ActivityHandle<?> handle = controller.addActivity(act);

        assertEquals(0, act.runCount);

        controller.start();
        controller.unpause();

        try {
            handle.future.get();
        } catch (InterruptedException | ExecutionException e) {
            fail(e);
        }

        assertEquals(20, act.runCount);
    }

    @Test
    void testActions() {
        Controller controller = new Controller();

        TestActivity act = new TestActivity();
        act.steps = 100;
        controller.setSleepPeriod(1);
        act.time = 1;
        ActivityHandle<?> handle = controller.addActivity(act);
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
        act.steps = 100;
        act.time = 1;

        ActivityHandle<?> handle = controller.addActivity(act);

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

    static class TestActivity extends Activity {
        volatile int runCount = 0;
        volatile int time = 0;
        private int steps;

        @Override
        public void restart() {

        }


        @Override
        public void runActivity(Controller.ControllerHandle handle) {
            for (int i = 0; i < this.steps; i++) {
                sleep(time);
                handle.yield(() -> runCount++);
            }
        }
    }

    static class TestAction implements IAction {
        volatile boolean ran = false;
        @Override
        public void actionPerformed() {
            ran = true;
        }
    }
}