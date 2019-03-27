import unittest
import itertools

import qiime2.core.type.poset as poset
from qiime2.core.type.tests.test_grammar import MockPredicate


class TestPOSet(unittest.TestCase):
    def test_sorting_complex(self):
        A = MockPredicate('A')
        B = MockPredicate('B')
        C = MockPredicate('C')
        D = MockPredicate('D')
        E = MockPredicate('E')

        for _, permutation in zip(range(250), itertools.permutations(
                [A, B, C & B & D, C & A & E, E & B, A & B & C & D, E,
                 B & A & E, C, D | E])):

            sort = list(poset.POSet(*permutation))
            ide = sort.index(D | E)
            ia = sort.index(A)
            ib = sort.index(B)
            ic = sort.index(C)
            ie = sort.index(E)
            ieb = sort.index(E & B)
            ibae = sort.index(B & A & E)
            icae = sort.index(C & A & E)
            icbd = sort.index(C & B & D)
            iabcd = sort.index(A & B & C & D)


            # verify topological sort
            for i in [ia, ic, ie, ide]:
                self.assertGreater(i, icae, sort)

            for i in [ia, ib, ie, ieb, ide]:
                self.assertGreater(i, ibae, sort)

            for i in [ib, ic, ide]:
                self.assertGreater(i, icbd, sort)

            for i in [ia, ib, ic, icbd, ide]:
                self.assertGreater(i, iabcd, sort)

            self.assertGreater(ide, ie, sort)


    def test_sorting_simple(self):
        A = MockPredicate('A')
        B = MockPredicate('B')
        C = MockPredicate('C')


        for permutation in itertools.permutations(
                [A, B, C, A & B, A & C, A & B & C]):

            sort = list(poset.POSet(*permutation))
            ia = sort.index(A)
            ib = sort.index(B)
            ic = sort.index(C)
            iab = sort.index(A & B)
            iac = sort.index(A & C)
            iabc = sort.index(A & B & C)

            for i in [ia, ib, ic, iab, iac]:
                self.assertGreater(i, iabc)

            for i in [ia, ic]:
                self.assertGreater(i, iac)

            for i in [ia, ib]:
                self.assertGreater(i, iab)


if __name__ == "__main__":
    unittest.main()
