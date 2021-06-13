import numpy as np
import matplotlib.pyplot as plt

days = np.array([20, 21, 23, 24, 25, 27, 29, 30], float)
death_counts = np.array([346, 362, 343, 339, 347, 346, 339, 394], float)

''''

*********************
  Lagrange Solution  
*********************

'''


def calculate_lagrange(predict, n, days, death_counts):
    f_x = 0.0
    for ii in range(n):
        prod = 1.0
        for jj in range(n):
            if jj == ii:
                continue
            else:
                prod *= (predict - days[jj]) / (days[ii] - days[jj])

        f_x += death_counts[ii] * prod
    return f_x


def draw_lagrange():
    x_cord_graph = np.linspace(days[0], days[-1])
    y_cord_graph = np.array([], float)

    for xp in x_cord_graph:
        yp = calculate_lagrange(xp, 8, days, death_counts)
        y_cord_graph = np.append(y_cord_graph, yp)

    plt.plot(days, death_counts, "ok", x_cord_graph, y_cord_graph)
    plt.show()


''''
*********************
  Direct Method Solution  
*********************
'''
def count_coefficients(x_val, count):
    coefficients = []
    for i in range(count):
        coefficients.append(x_val ** i)
    return coefficients

def create_coefficients_matrix(points):
    new_matrix = np.zeros((len(points), len(points)))
    count = 0
    for x_val in points:
        new_matrix[count] = count_coefficients(x_val, len(points))
        count += 1
    return new_matrix

def calculate_direct_method(x, y, n, xn):
    x = np.array(x)
    y = np.array(y).reshape(n, 1)
    yn = 0

    c_matrix = create_coefficients_matrix(x)
    s_matrix = np.linalg.solve(c_matrix, y)

    for i in range(len(s_matrix)):
        yn += (s_matrix[i]) * (pow(xn, i))

    x = np.append(x, xn)
    x.sort()

    index_xn = list(x).index(26)
    y = np.insert(y, index_xn, yn)

    plt.plot(days, death_counts, "ok", x, y)
    plt.show()

    return yn


''''

*********************
  Newton’s Divided Difference Solution  
*********************

'''


def coefficient(x_array, y_array):
    x_array.astype(float)
    y_array.astype(float)
    n = len(x_array)

    coefficient_array = []

    for i in range(n):
        coefficient_array.append(y_array[i])

    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coefficient_array[i] = float(coefficient_array[i] - coefficient_array[i - 1]) / float(
                x_array[i] - x_array[i - j])

    return np.array(coefficient_array)


def evaluate(coefficient_array, x_array, p):
    x_array.astype(float)
    n = len(coefficient_array) - 1
    fx = coefficient_array[n]

    for i in range(n - 1, -1, -1):
        fx = fx * (p - x_array[i]) + coefficient_array[i]

    return fx


def calculate_NDDD(days, deaths, p):
    coff = coefficient(days, deaths)

    result_pol = np.polynomial.Polynomial([0.])
    coffnumber = coff.shape[0]

    for i in range(coffnumber):
        p = np.polynomial.Polynomial([1.])
        for j in range(i):
            p_temp = np.polynomial.Polynomial([-days[j], 1.])
            p = np.polymul(p, p_temp)
        p *= coff[i]
        result_pol = np.polyadd(result_pol, p)

    coefficients = np.flip(result_pol[0].coef)

    y_axis = np.polyval(coefficients, days)

    plt.plot(days, y_axis, "g")
    plt.show()

    return evaluate(coff, days, p)


'''
 TEST / PRINT 
'''

'''
 # Lagrange
'''
# lagrange = calculate_lagrange(26, 8, days, death_counts)
# print("Result is " + str(lagrange))
# draw_lagrange()


'''
 # Direct Method of Polynomial Interpolation
'''
# directMethod = calculate_direct_method(days, death_counts, 8, 26)
# print("Result is " + str(directMethod))


'''
 # Newton’s Divided Difference Interpolation
'''
# result = calculate_NDDD(days, death_counts, 26)
# print("Result is " + str(result))
