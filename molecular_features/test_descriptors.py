import unittest
import descriptors
import filecmp


class TestCalc(unittest.TestCase):
    def test_compute_descriptors(self):
        descriptors.compute_descriptors("data/output1.json", "data/result1.csv", False)
        same = filecmp.cmp("data/true1.csv", "data/result1.csv")
        self.assertEqual(same, True)

        descriptors.compute_descriptors("data/output2.json", "data/result2.csv", True)
        same = filecmp.cmp("data/true2.csv", "data/result2.csv")
        self.assertEqual(same, True)

        descriptors.compute_descriptors("data/output3.json", "data/result3.csv", False)
        same = filecmp.cmp("data/true3.csv", "data/result3.csv")
        self.assertEqual(same, True)


if __name__ == "__main__":
    unittest.main()
