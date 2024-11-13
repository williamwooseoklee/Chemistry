import numpy as np
import string
import sympy.matrices as sp

def multiplyElements(compounds):
    new_compounds = []
    for compound in compounds:
        added_compound = compound
        added_string = ""
        added = False
        for i in range(len(compound)-1):
            if compound[i + 1].isnumeric():
                added_string += (int(compound[i + 1]) - 1) * compound[i]
                added_compound = compound + added_string
        new_compounds.append(added_compound)

    return new_compounds

def countElements(compounds, elements):

    totalelementcount = []
    for new_compound in compounds:
        element_count = []
        for element in elements:
            element_count.append(new_compound.count(element))
        totalelementcount.append(element_count)
    totalelementcount = np.array(totalelementcount)
    return np.transpose(totalelementcount)

def matrix_augmentation(coefficient_matrix):
    rows = coefficient_matrix.shape[0]
    columns = coefficient_matrix.shape[1]

    if rows == columns:
        coefficient_matrix = np.float64(sp.Matrix(coefficient_matrix).echelon_form())

    rank = np.linalg.matrix_rank(coefficient_matrix)
    nullity = columns - rank

    if nullity != 0:
        if rows == columns:
            coefficient_matrix = coefficient_matrix[0:(rows - nullity), :]
        for row_n in range(nullity):
            new_row = np.zeros(columns)
            new_row[-1 * (nullity - row_n)] = 1
            coefficient_matrix = np.vstack((coefficient_matrix, new_row))

    return [coefficient_matrix,nullity]

def matrix_to_equation(coefficient_matrix, nullity):
    if nullity != 0:
        coefficient_matrix_inv = np.linalg.inv(coefficient_matrix)
        solution_coefficients = coefficient_matrix_inv[:, -1]
        solution_coefficients = solution_coefficients / min(abs(np.ma.masked_equal(solution_coefficients, 0.0, copy=False)))

        output = ""
        species = np.array(reactants + products)
        equal_used = False
        for i in range(len(solution_coefficients)):
            output += (str(round(abs(solution_coefficients[i]))) + str(species[i]))

            if i < len(solution_coefficients) - 1:
                if solution_coefficients[i + 1] > 0 and not equal_used:
                    output += " = "
                    equal_used = True
                else:
                    output += " + "
    else:
        return "Equation cannot be balanced."

    return output

input = input("Enter equation: ")

elements = []

input = input.replace(" ","")

currentElement = ""

for i in range(len(input)-1):
    if input[i].isalpha() and (input[i+1].isdigit() or input[i+1].isupper() or (input[i+1] in string.punctuation)):
        currentElement += input[i]
        if currentElement in elements:
            currentElement = ""
        else:
            elements.append(currentElement)
            currentElement = ""
    elif input[i].isalpha():
        currentElement += input[i]

input = input.split("=")
reactants = input[0].split("+")
products = input[1].split("+")

coefficient_matrix = np.hstack((countElements(multiplyElements(reactants),elements),countElements(multiplyElements(products),elements)))

augmentation = matrix_augmentation(coefficient_matrix)
coefficient_matrix = augmentation[0]
nullity = augmentation[1]

output = matrix_to_equation(coefficient_matrix,nullity)



print(output)