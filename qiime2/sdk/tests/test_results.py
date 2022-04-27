# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pickle
import unittest

from qiime2.sdk import Results


class TestResults(unittest.TestCase):
    def test_tuple_subclass(self):
        self.assertTrue(issubclass(Results, tuple))
        self.assertIsInstance(Results(['a', 'b'], [42, 43]), tuple)

    def test_tuple_cast(self):
        r = Results(['a', 'b'], [42, 43])

        t = tuple(r)

        self.assertIs(type(t), tuple)
        self.assertEqual(t, (42, 43))

    def test_callable_return_and_unpacking(self):
        def f():
            return Results(['a', 'b'], [42, 43])

        a, b = f()

        self.assertEqual(a, 42)
        self.assertEqual(b, 43)

    def test_constructor_iterable(self):
        r = Results(iter(['a', 'b']), iter([42, 43]))

        self.assertEqual(tuple(r), (42, 43))
        self.assertEqual(r._fields, ('a', 'b'))

    def test_constructor_len_mismatch(self):
        with self.assertRaises(ValueError):
            Results(['a', 'b'], [42])

    def test_pickle(self):
        r = Results(['a', 'b'], [42, 'abc'])

        pickled = pickle.dumps(r)
        unpickled = pickle.loads(pickled)

        self.assertEqual(unpickled, r)

    def test_field_attributes(self):
        r = Results(['foo', 'bar'], [42, 'abc'])

        self.assertEqual(r.foo, 42)
        self.assertEqual(r.bar, 'abc')

        with self.assertRaises(AttributeError):
            r.baz

    def test_per_instance_field_attributes(self):
        # Field attributes are added to a `Results` instance, not the type.
        r1 = Results(['foo', 'bar'], [42, 'abc'])
        r2 = Results(['x'], [42.0])

        for attr in 'foo', 'bar', 'x':
            self.assertFalse(hasattr(Results, attr))

        self.assertTrue(hasattr(r1, 'foo'))
        self.assertTrue(hasattr(r1, 'bar'))
        self.assertFalse(hasattr(r1, 'x'))

        self.assertFalse(hasattr(r2, 'foo'))
        self.assertFalse(hasattr(r2, 'bar'))
        self.assertTrue(hasattr(r2, 'x'))

    def test_index_access(self):
        r = Results(['foo', 'bar'], [42, 'abc'])

        self.assertEqual(r[0], 42)
        self.assertEqual(r[1], 'abc')

        with self.assertRaises(IndexError):
            r[2]

    def test_immutability(self):
        r = Results(['foo', 'bar'], [42, 'abc'])

        # Setter for existing attribute.
        with self.assertRaises(AttributeError):
            r.bar = 999

        # Setter for new attribute.
        with self.assertRaises(AttributeError):
            r.baz = 999

        # Deleter for existing attribute.
        with self.assertRaises(AttributeError):
            del r.bar

        # Deleter for new attribute.
        with self.assertRaises(AttributeError):
            del r.baz

        with self.assertRaises(TypeError):
            r[0] = 999

    def test_eq_same_obj(self):
        r = Results(['a', 'b'], [1, 2])

        self.assertEqual(r, r)

    def test_eq_subclass(self):
        class ResultsSubclass(Results):
            pass

        r1 = Results(['foo'], ['abc'])
        r2 = ResultsSubclass(['foo'], ['abc'])

        self.assertEqual(r1, r2)

    def test_eq_different_source_types(self):
        r1 = Results(iter(['a', 'b']), iter([42, 43]))
        r2 = Results(['a', 'b'], [42, 43])

        self.assertEqual(r1, r2)

    def test_eq_empty(self):
        r1 = Results([], [])
        r2 = Results([], [])

        self.assertEqual(r1, r2)

    def test_eq_nonempty(self):
        r1 = Results(['foo', 'bar'], ['abc', 'def'])
        r2 = Results(['foo', 'bar'], ['abc', 'def'])

        self.assertEqual(r1, r2)

    def test_ne_type(self):
        r1 = Results(['foo', 'bar'], ['abc', 'def'])
        r2 = ('abc', 'def')

        self.assertNotEqual(r1, r2)

    def test_ne_fields(self):
        r1 = Results(['foo', 'bar'], ['abc', 'def'])
        r2 = Results(['foo', 'baz'], ['abc', 'def'])

        self.assertNotEqual(r1, r2)

    def test_ne_values(self):
        r1 = Results(['foo', 'bar'], ['abc', 'def'])
        r2 = Results(['foo', 'bar'], ['abc', 'xyz'])

        self.assertNotEqual(r1, r2)

    def test_repr_empty(self):
        r = Results([], [])

        self.assertTrue(repr(r).startswith('Results'))
        self.assertTrue(repr(r).endswith('---'))

    def test_repr_single(self):
        r = Results(['a'], [42])

        self.assertTrue(repr(r).startswith('Results'))
        self.assertTrue(repr(r).endswith('a = 42'))

    def test_repr_multiple(self):
        r = Results(['a', 'foo'], [42, 'abc'])

        self.assertTrue(repr(r).startswith('Results'))
        self.assertTrue('a   = 42' in repr(r))
        self.assertTrue(repr(r).endswith("foo = 'abc'"))


if __name__ == '__main__':
    unittest.main()
