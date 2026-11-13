import argparse
from fractions import Fraction
import json
from math import comb, gcd, lcm
from pathlib import Path


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def polynomial_add(first, second, second_scale=1):
    result = [0] * max(len(first), len(second))
    for degree, coefficient in enumerate(first):
        result[degree] += coefficient
    for degree, coefficient in enumerate(second):
        result[degree] += second_scale * coefficient
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def polynomial_sequence(initial, following, maximum):
    values = [initial, following]
    for index in range(2, maximum + 1):
        values.append(polynomial_add(values[index - 1], [0] + values[index - 2], -1))
    return values


def expected_row(index):
    quotient, residue = divmod(index, 6)
    values = [Fraction(0), Fraction(1, 3 * quotient + 1),
              Fraction(1, 2 * quotient + 1), Fraction(1, 3 * quotient + 2),
              Fraction(0), Fraction(-1, 6 * quotient + 6)]
    return values[residue]


def matrix_product(first, second, modulus):
    return tuple(sum(first[2 * row + offset] * second[2 * offset + column]
                     for offset in range(2)) % modulus
                 for row in range(2) for column in range(2))


def check():
    maximum = 240
    triangle = [[Fraction(0)] * (maximum - row + 1) for row in range(maximum + 1)]
    for row in range(maximum):
        triangle[row][1] = Fraction((-1) ** row, row + 1)
    for distance in range(2, maximum + 1):
        for row in range(maximum - distance + 1):
            triangle[row][distance] = triangle[row][distance - 1] + triangle[row + 1][distance - 2]
    values = triangle[0]
    denominators = [value.denominator for value in values]
    numerator_block = [0, 1, 1, 1, 0, -1]
    numerator_coefficients = [1, 1, 1, 2, 1, 6, -1, 2, 1, 1, -1]
    for index, value in enumerate(values):
        require(value == expected_row(index), ('first_row', index))
        quotient, residue = divmod(index, 6)
        expected_denominator = [1, 3 * quotient + 1, 2 * quotient + 1,
                                3 * quotient + 2, 1, 6 * quotient + 6][residue]
        require(value.denominator == expected_denominator, ('denominator', index))
        require(value.numerator == numerator_block[residue], ('numerator', index))
        logarithmic = Fraction(0)
        if (index + 1) % 2 == 0:
            logarithmic += Fraction(2, index + 1)
        if (index + 1) % 3 == 0:
            logarithmic += Fraction(3 * (-1) ** ((index + 1) // 3 + 1), index + 1)
        require(value == logarithmic, ('logarithmic_coefficient', index))
        coefficient = denominators[index]
        if index >= 6:
            coefficient -= 2 * denominators[index - 6]
        if index >= 12:
            coefficient += denominators[index - 12]
            require(denominators[index] == 2 * denominators[index - 6] - denominators[index - 12],
                    ('denominator_recurrence', index))
        expected = numerator_coefficients[index] if index < len(numerator_coefficients) else 0
        require(coefficient == expected, ('denominator_generating_function', index))

    numerators = [value.numerator for value in values]
    quotients = [Fraction(index + 1, denominator) for index, denominator in enumerate(denominators)]
    unsigned_triangle = [[Fraction(0)] * (maximum - row + 1) for row in range(maximum + 1)]
    for row in range(maximum):
        unsigned_triangle[row][1] = Fraction(1, row + 1)
    for distance in range(2, maximum + 1):
        for row in range(maximum - distance + 1):
            unsigned_triangle[row][distance] = (unsigned_triangle[row][distance - 1]
                                                + unsigned_triangle[row + 1][distance - 2])
    lucas = [2, 1]
    for index in range(2, 707):
        lucas.append(lucas[-1] + lucas[-2])
    for index in range(maximum + 1):
        require(unsigned_triangle[0][index] == Fraction(lucas[index + 1] - 1, index + 1),
                ('unsigned_lucas_formula', index))
    for index in range(120):
        require(sum(unsigned_triangle[row][2 * index + 1 - 2 * row] for row in range(index + 1))
                == unsigned_triangle[0][2 * index + 2], ('unsigned_antidiagonal', index))
    for index in range(41):
        transformed = sum(((-1) ** (index - term) * comb(index, term) * unsigned_triangle[0][term]
                           for term in range(index + 1)), Fraction(0))
        require(transformed == (-1) ** (index + 1) * unsigned_triangle[0][index],
                ('unsigned_binomial_transform', index))
    for index in range(81):
        half = (index + 1) // 2
        expected = lcm(*range(1, half + 1))
        for label, array in [('alternating', triangle), ('unsigned', unsigned_triangle)]:
            actual = lcm(*(array[row][index - 2 * row].denominator for row in range(half)))
            require(actual == expected, ('antidiagonal_denominator_lcm', label, index))
    modulus = 705
    powers = {1: (1, 1, 1, 0)}
    for index in range(1, 10):
        powers[2 ** index] = matrix_product(powers[2 ** (index - 1)], powers[2 ** (index - 1)], modulus)
    require(powers[64] == (565, 423, 423, 142), 'matrix_64')
    require(powers[128] == (424, 141, 141, 283), 'matrix_128')
    require(powers[512] == (283, 564, 564, 424), 'matrix_512')
    product = (1, 0, 0, 1)
    for exponent in [1, 64, 128, 512]:
        product = matrix_product(product, powers[exponent], modulus)
    require(product == (283, 565, 565, 423) and lucas[705] % modulus == 1,
            'composite_integer_rank')
    require(modulus == 3 * 5 * 47 and (lucas[705] - 1) // modulus > 0,
            'composite_rank_factorization_and_positivity')
    quotient_coefficients = [1, 2, 3, 2, 5, 1, 5, -2, -3, -2, 1, -1]
    for index in range(maximum + 1):
        quotient, residue = divmod(index, 6)
        require(quotients[index] == [6 * quotient + 1, 2, 3, 2, 6 * quotient + 5, 1][residue],
                ('quotient_tuple', index))
        require(quotients[index].denominator == 1, ('quotient_integrality', index))
        coefficient = quotients[index]
        if index >= 6:
            coefficient -= 2 * quotients[index - 6]
        if index >= 12:
            coefficient += quotients[index - 12]
        expected = quotient_coefficients[index] if index < len(quotient_coefficients) else 0
        require(coefficient == expected, ('quotient_generating_function', index))
        sign = -1 if ((index - 1) // 3) % 2 else 1
        floor_value = (1 + sign) // 2 - (index + 1) // 6 + index // 6
        require(numerators[index] == floor_value, ('numerator_floor', index))
        coefficient = sum((-1) ** shift * numerators[index - shift]
                          for shift in range(min(5, index) + 1))
        require(coefficient == {1: 1, 3: 1, 4: -1}.get(index, 0),
                ('numerator_generating_function_and_recurrence', index))

    def auxiliary(entry, index):
        if entry == 'A130196':
            return 1 if index % 3 == 0 else 2
        if entry == 'A052901':
            return 3 if index % 3 == 0 else 2
        if entry == 'A047245':
            return 6 * ((index - 1) // 3) + (index - 1) % 3 + 1
        if entry == 'A047238':
            return 6 * ((index - 1) // 2) + (2 if index % 2 == 0 else 0)
        if entry == 'A026741':
            return index if index % 2 else index // 2
        if entry == 'A045896':
            return Fraction(index, (index + 1) * (index + 2)).denominator
        if entry == 'A051176':
            return index // gcd(index, 3)
        if entry == 'A146535':
            return Fraction(2 * index - 1, 3).numerator
        if entry == 'A306368':
            return Fraction((index + 3) * (index + 4), (index + 1) * (index + 2)).numerator
        if entry == 'A194767':
            return Fraction(2 * index + 1, (index + 1) * (index + 2)).denominator
        if entry == 'A139550':
            return lcm(*range(1, index // 2 + 1))
        if entry == 'A189731':
            return Fraction(lucas[index + 1] - 1, index + 1).numerator
        if entry == 'A242926':
            return Fraction(lucas[index + 1] - 1, index + 1).denominator
        if entry in {'A175385', 'A175386'}:
            fraction = sum((Fraction(comb(2 * index - term - 1, term - 1), term)
                            for term in range(1, index + 1)), Fraction(0))
            return fraction.numerator if entry == 'A175385' else fraction.denominator
        raise ValueError(entry)

    for index in range(80):
        require(quotients[2 * index + 1] == auxiliary('A130196', index + 1),
                ('odd_subsequence', index))
        require(quotients[auxiliary('A047245', index + 1)] == auxiliary('A052901', index + 2),
                ('ranked_subsequence', index))
        require(sum(quotients[3 * index:3 * index + 3]) == auxiliary('A047238', index + 3),
                ('triple_sums', index))

    source_data = json.loads(Path(__file__).with_name('source-terms.json').read_text())
    source_comparisons = {}
    for entry, record in source_data.items():
        for position, expected in enumerate(record['terms']):
            index = position + record['offset']
            if entry == 'A189733':
                actual = denominators[index]
            elif entry == 'A278718':
                actual = numerators[index]
            elif entry == 'A267068':
                actual = quotients[index]
            else:
                actual = auxiliary(entry, index)
            require(actual == expected, ('source_terms', entry, index))
        source_comparisons[entry] = len(record['terms'])
    for index in range(80):
        require(unsigned_triangle[0][2 * index + 1].numerator == auxiliary('A175385', index + 1),
                ('unsigned_odd_numerator', index))
        require(unsigned_triangle[0][2 * index + 1].denominator == auxiliary('A175386', index + 1),
                ('unsigned_odd_denominator', index))
    require(unsigned_triangle[0][3].denominator == 2 and auxiliary('A175386', 1) == 1,
            'unshifted_unsigned_denominator_counterexample')

    second_denominators = [auxiliary('A194767', index) for index in range(maximum + 1)]
    second_coefficients = [2, 2, 12, 14, 4, 6, 2]
    for index, denominator in enumerate(second_denominators):
        quotient, residue = divmod(index, 3)
        expected = [(3 * quotient + 1) * (3 * quotient + 2),
                    (quotient + 1) * (3 * quotient + 2),
                    3 * (quotient + 1) * (3 * quotient + 4)][residue]
        require(denominator == expected, ('second_kind_residues', index))
        require(denominator == (index + 1) * auxiliary('A051176', index + 2),
                ('second_kind_first_alias', index))
        require(auxiliary('A194767', index + 2) == 2 * auxiliary('A306368', index),
                ('second_kind_second_alias', index))
        expected = second_coefficients[index] if index < len(second_coefficients) else 0
        actual = sum(scale * second_denominators[index - shift]
                     for shift, scale in [(0, 1), (3, -3), (6, 3), (9, -1)] if index >= shift)
        require(actual == expected, ('second_kind_generating_function', index))
    for index in range(81):
        diagonal = 2 * triangle[index][1] - triangle[index][0]
        first = 2 * triangle[index][2] - triangle[index][1]
        second = 2 * triangle[index][3] - triangle[index][2]
        third = 2 * triangle[index][4] - triangle[index][3]
        require(diagonal == Fraction(2 if index % 2 == 0 else -1,
                                     auxiliary('A026741', index + 1)), ('second_kind_main', index))
        require(first == Fraction((-1) ** index, index + 1), ('second_kind_first', index))
        require(second == Fraction((-1) ** (index + 1) * auxiliary('A026741', index),
                                   auxiliary('A045896', index)), ('second_kind_second', index))
        require(third.denominator == second_denominators[index], ('second_kind_denominator', index))
        require(third.numerator == (-1) ** (index + 1) * auxiliary('A146535', index + 1),
                ('corrected_second_kind_numerator', index))
    require(auxiliary('A146535', 2) == 1 and 2 * triangle[2][4] - triangle[2][3] == Fraction(-5, 12),
            'unshifted_numerator_counterexample')

    square_maximum = 40
    square = [[Fraction(0)] * (square_maximum + 1) for row in range(square_maximum + 1)]
    for row in range(square_maximum + 1):
        for column in range(row, square_maximum + 1):
            square[row][column] = triangle[row][column - row]
        for column in range(row):
            square[row][column] = square[row - 1][column + 1] - square[row - 1][column]
    for row in range(square_maximum + 1):
        for column in range(square_maximum + 1):
            require(square[row][column] == (-1) ** (row + column + 1) * square[column][row],
                    ('reflection', row, column))
    for row in range(square_maximum):
        for column in range(square_maximum):
            require(square[row + 1][column] == square[row][column + 1] - square[row][column],
                    ('difference_identity', row, column))
    for index in range(square_maximum + 1):
        transformed = sum(((-1) ** (index - term) * comb(index, term) * values[term]
                           for term in range(index + 1)), Fraction(0))
        require(square[index][0] == transformed == (-1) ** (index + 1) * values[index],
                ('binomial_transform', index))

    polynomial_maximum = 80
    primary = polynomial_sequence([0], [1], polynomial_maximum + 1)
    companion = polynomial_sequence([2], [1], polynomial_maximum + 1)
    for index in range(polynomial_maximum + 1):
        derivative = [(degree + 1) * coefficient
                      for degree, coefficient in enumerate(companion[index + 1][1:])]
        require(polynomial_add(derivative or [0], primary[index], index + 1) == [0],
                ('companion_derivative', index))
        for row in range(11):
            moment = sum((Fraction(coefficient, row + degree + 1)
                          for degree, coefficient in enumerate(primary[index])), Fraction(0))
            require(triangle[row][index] == (-1) ** row * moment,
                    ('polynomial_moment', row, index))
            if index >= 1:
                explicit = (-1) ** row * sum(
                    (Fraction((-1) ** degree * comb(index - 1 - degree, degree), row + degree + 1)
                     for degree in range((index - 1) // 2 + 1)), Fraction(0))
                require(triangle[row][index] == explicit, ('all_entry_sum', row, index))
        endpoint = Fraction(1 - sum(companion[index + 1]), index + 1)
        require(values[index] == endpoint, ('companion_endpoint', index))

    return {
        'status': 'pass',
        'arithmetic': 'exact integers and fractions.Fraction',
        'first_row': {'indices': [0, maximum], 'comparisons': maximum + 1},
        'reduced_denominators': {'indices': [0, maximum], 'comparisons': maximum + 1},
        'reduced_numerators': {'indices': [0, maximum], 'period_block': numerator_block},
        'logarithmic_coefficients': {'indices': [0, maximum]},
        'denominator_recurrence': {'indices': [12, maximum], 'comparisons': maximum - 11},
        'denominator_generating_function': {'coefficients': [0, maximum]},
        'numerator_floor_and_generating_function': {'indices': [0, maximum]},
        'numerator_order_five_recurrence': {'indices': [5, maximum]},
        'quotient_tuple_integrality_and_generating_function': {'indices': [0, maximum]},
        'quotient_recurrence': {'indices': [12, maximum]},
        'quotient_sequence_identifications': {'indices': [0, 79], 'comparisons': 240},
        'published_source_terms': source_comparisons,
        'unsigned_lucas_formula': {'indices': [0, maximum]},
        'unsigned_binomial_transform': {'indices': [0, 40]},
        'unsigned_antidiagonal': {'indices': [0, 119]},
        'antidiagonal_denominator_lcm': {'antidiagonals': [0, 80], 'arrays': ['alternating', 'unsigned'],
                                        'comparisons': 162},
        'unsigned_odd_bisections': {'indices': [0, 79], 'comparisons': 160,
                                   'corrected_denominator_counterexample': {'index': 3,
                                                                          'actual': 2, 'unshifted': 1}},
        'composite_integer_rank': {'rank': modulus, 'zero_based_index': 704, 'factorization': [3, 5, 47],
                                   'lucas_remainder': lucas[705] % modulus,
                                   'matrix_power_residues': {str(exponent): list(powers[exponent])
                                                             for exponent in [64, 128, 512]},
                                   'matrix_705_residue': list(product)},
        'second_kind_residues_generating_function_and_aliases': {'indices': [0, maximum]},
        'second_kind_first_four_diagonals': {'indices': [0, 80], 'comparisons': 324},
        'corrected_second_kind_numerator': {'indices': [0, 80],
                                           'counterexample_to_unshifted_formula': {'index': 2, 'value': '-5/12',
                                                                                  'unshifted_numerator': 1}},
        'square_reflection': {'rows': [0, square_maximum], 'columns': [0, square_maximum],
                              'comparisons': (square_maximum + 1) ** 2},
        'square_difference_identity': {'rows': [0, square_maximum - 1],
                                      'columns': [0, square_maximum - 1]},
        'binomial_transform': {'indices': [0, square_maximum], 'comparisons': square_maximum + 1},
        'polynomial_moments': {'rows': [0, 10], 'indices': [0, polynomial_maximum],
                               'comparisons': 11 * (polynomial_maximum + 1)},
        'companion_derivative_and_endpoint': {'indices': [0, polynomial_maximum]},
        'all_entry_binomial_sum': {'rows': [0, 10], 'distances': [1, polynomial_maximum],
                                  'comparisons': 11 * polynomial_maximum},
        'scope': 'Finite exact checks of the stated ranges; universal proofs are in note.tex.'
    }


def main():
    parser = argparse.ArgumentParser(description='Exact finite checks for the OEIS A189733 formulas.')
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    encoded = json.dumps(check(), indent=2, sort_keys=True) + '\n'
    if arguments.output:
        arguments.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
