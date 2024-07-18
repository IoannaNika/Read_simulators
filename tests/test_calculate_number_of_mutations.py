import unittest
from src.calculate_number_of_mutations import calc_n_mutations

class TestCalculateNumberOfMutations(unittest.TestCase):
    def test_calc_n_mutations1(self):
        read1 = "AC"
        read2 = "A"
        self.assertEqual(calc_n_mutations(read1, read2), 0)
    
    def test_calc_n_mutations2(self):
        read1 = "AC"
        read2 = "AC"
        self.assertEqual(calc_n_mutations(read1, read2), 0)

    def test_calc_n_mutations3(self):
        read1 = "AC"
        read2 = "AT"
        self.assertEqual(calc_n_mutations(read1, read2), 1)

    def test_calc_n_mutations4(self):
        read1 = "C"
        read2 = "AG"
        self.assertEqual(calc_n_mutations(read1, read2), 1)

if __name__ == '__main__':
    unittest.main()