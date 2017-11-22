package etomica.util.collections;

import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.*;

public class IndexMapTest {
    public IndexMap<String> map;

    @Before
    public void setUp() {
        this.map = new IndexMap<>();
    }

    @Test
    public void testSize() {
        assertEquals(0, map.size());

        map.put(0, "test");
        assertEquals(1, map.size());

        map.put(100, "test2");
        assertEquals(2, map.size());
    }

    @Test
    public void testIsEmpty() {
        assertTrue(map.isEmpty());
        map.put(1, "test");
        assertFalse(map.isEmpty());
        map.remove(1);
        assertTrue(map.isEmpty());
    }

    @Test
    public void containsKey() {
        assertFalse(map.containsKey(0));

        map.put(0, "");
        map.put(5, "");
        assertTrue(map.containsKey(5));
        assertFalse(map.containsKey(3));
    }

    @Test
    public void containsValue() {
        assertFalse(map.containsValue("test"));
        map.put(0, "foo");
        map.put(20, "test");
        assertTrue(map.containsValue("test"));
    }

    @Test
    public void testGet() {
        assertNull(map.get(0));

        for(int i = 0; i < 100; i++) {
            map.put(i, String.valueOf(i));
        }

        map.put(105, "105");

        assertEquals("50", map.get(50));
        assertEquals("105", map.get(105));
        assertNull(map.get(104));
        assertNull(map.get(-1));
        assertNull(map.get(106));
    }

    @Test(expected = IllegalArgumentException.class)
    public void testPutNullValue() {
        map.put(5, null);
    }

    @Test
    public void testRemove() {
        assertNull(map.remove(0));
        assertEquals(0, map.size());

        map.put(0, "0");
        map.put(1, "1");
        map.put(5, "5");
        assertEquals("1", map.remove(1));
        assertEquals(2, map.size());

        map.put(90, "90");
        assertEquals("5", map.remove(5));
        assertEquals(2, map.size());
    }

    @Test
    public void testPutAll() {
        Map<Integer, String> map2 = new HashMap<>();
        for(int i = 0; i < 100; i++) {
            map2.put(i, String.valueOf(i));
        }

        map.putAll(map2);
        assertEquals(100, map.size());
        assertEquals("40", map.get(40));
    }

    @Test
    public void testClear() {
        map.put(50, "50");
        map.clear();
        assertEquals(0, map.size());
    }

    @Test(expected = UnsupportedOperationException.class)
    public void testKeySet() {
        map.keySet();
    }

    @Test()
    public void testValues() {
        Collection<String> c = map.values();
        assertTrue(c.isEmpty());

        map.put(0, "0");
        map.put(1, "1");
        map.put(5, "5");
        assertTrue(c.containsAll(Arrays.asList("0", "1", "5")));

        map.remove(1);
        assertTrue(c.containsAll(Arrays.asList("0", "5")));
        assertFalse(c.contains("1"));

        // test resizing array
        map.put(100, "100");
        assertTrue(c.containsAll(Arrays.asList("0", "5", "100")));
    }

    @Test(expected = UnsupportedOperationException.class)
    public void testEntrySet() {
        map.entrySet();
    }
}
