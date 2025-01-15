import numpy as np


def fit_curve(points,n=10,progress_callback=None):
    "n correspond au nombre de points par sous ensemble"
    if len(points) < 2:
        return []
    beziers = []
    for i in range(0, len(points) - 1, n-1):
        subset_points = points[i:i + n]
        if len(subset_points) < 2:
            break

        # Calcul des tangentes aux extr�mit�s (vecteurs unitaires)
        left_tangent = create_tangent(subset_points[1], subset_points[0])
        right_tangent = create_tangent(subset_points[-2], subset_points[-1])
        #calcule les béziers pour chaque sous ensemble
        beziers += fit_cubic(subset_points, left_tangent, right_tangent, progress_callback)

       
    return beziers
# Fonction pour ajuster une courbe de B�zier cubique � un ensemble de points
def fit_cubic(points, left_tangent, right_tangent, progress_callback=None):
    max_iterations = 1000

    # Param�trisation initiale bas�e sur la longueur des cordes
    u = chord_length_parameterize(points)
    bez_curve = generate(points, u, u, left_tangent, right_tangent, progress_callback)

    #
    u_prime = u
    for i in range(max_iterations):
            u_prime = reparameterize(bez_curve, points, u_prime)
            bez_curve = generate_bezier(points, params_prime, left_tangent, right_tangent)


    return [bez_curve]

# G�n�rer une courbe de B�zier cubique
def generate_bezier(points, parameters, left_tangent, right_tangent):
    first_point = points[0]
    last_point = points[-1]
    bez_curve = [first_point, None, None, last_point]
    

    # Calculer les coefficients de la matrice A
    A = np.zeros((len(parameters), 2, 2))
    for i, u in enumerate(parameters):
        ux = 1 - u
        A[i][0] = left_tangent * (3 * u * (ux ** 2))
        A[i][1] = right_tangent * (3 * ux * (u ** 2))

    # R�soudre le syst�me lin�aire pour trouver les points de contr�le interm�diaires
    C = np.zeros((2, 2))
    X = np.zeros(2)
    for i, u in enumerate(parameters):
        a = A[i]
        C[0][0] += np.dot(a[0], a[0])
        C[0][1] += np.dot(a[0], a[1])
        C[1][0] += np.dot(a[0], a[1])
        C[1][1] += np.dot(a[1], a[1])
        tmp = np.subtract(points[i], bezier_q([first_point, first_point, last_point, last_point], u))
        X[0] += np.dot(a[0], tmp)
        X[1] += np.dot(a[1], tmp)

    # Calcul des d�terminants pour obtenir les valeurs alpha
    det_C0_C1 = (C[0][0] * C[1][1]) - (C[1][0] * C[0][1])
    det_C0_X = (C[0][0] * X[1]) - (C[1][0] * X[0])
    det_X_C1 = (X[0] * C[1][1]) - (X[1] * C[0][1])

    alpha_l = 0 if det_C0_C1 == 0 else det_X_C1 / det_C0_C1
    alpha_r = 0 if det_C0_C1 == 0 else det_C0_X / det_C0_C1
 
        
    # Si alpha est trop petit, utiliser une approximation simple
    seg_length = np.linalg.norm(np.subtract(first_point, last_point))
    alpha_l = max(0, min(alpha_l, seg_length / 1.5))  # Limite à la moitié de la longueur du segment
    alpha_r = max(0, min(alpha_r, seg_length / 1.5))
    epsilon = 1.0e-6 * seg_length
    if alpha_l < epsilon or alpha_r < epsilon:
        bez_curve[1] = np.add(first_point, left_tangent * (seg_length / 3.0))
        bez_curve[2] = np.add(last_point, right_tangent * (seg_length / 3.0))
    else:
        bez_curve[1] = np.add(first_point, left_tangent * alpha_l)
        bez_curve[2] = np.add(last_point, right_tangent * alpha_r)

    return bez_curve

# R�ajuster la param�trisation des points pour am�liorer l'ajustement
def reparameterize(bezier, points, parameters):
    return [newton_raphson_root_find(bezier, points[i], u) for i, u in enumerate(parameters)]

# Utiliser la m�thode de Newton-Raphson pour trouver une meilleure param�trisation
def newton_raphson_root_find(bez, point, u):
    d = np.subtract(bezier_q(bez, u), point)
    qprime = bezier_qprime(bez, u)
    numerator = np.dot(d, qprime)
    denominator = np.sum(np.square(qprime)) + 2 * np.dot(d, bezier_qprimeprime(bez, u))

    if denominator == 0:
        return u
    else:
        return u - (numerator / denominator)

# Param�trisation initiale bas�e sur la longueur des cordes
def chord_length_parameterize(points):
    u = [0]
    for i in range(1, len(points)):
        u.append(u[-1] + np.linalg.norm(np.subtract(points[i], points[i - 1])))
    return [x / u[-1] for x in u]



# �valuer la courbe de B�zier au param�tre t
def bezier_q(ctrl_poly, t):
    tx = 1.0 - t
    pA = np.multiply(ctrl_poly[0], (tx ** 3))
    pB = np.multiply(ctrl_poly[1], (3 * (tx ** 2) * t))
    pC = np.multiply(ctrl_poly[2], (3 * tx * (t ** 2)))
    pD = np.multiply(ctrl_poly[3], (t ** 3))
    return pA + pB + pC + pD

# Calculer la premi�re d�riv�e de la courbe de B�zier
def bezier_qprime(ctrl_poly, t):
    tx = 1.0 - t
    pA = (ctrl_poly[1] - ctrl_poly[0]) * (3 * (tx ** 2))
    pB = (ctrl_poly[2] - ctrl_poly[1]) * (6 * tx * t)
    pC = (ctrl_poly[3] - ctrl_poly[2]) * (3 * (t ** 2))
    return pA + pB + pC

# Calculer la seconde d�riv�e de la courbe de B�zier
def bezier_qprimeprime(ctrl_poly, t):
    return (ctrl_poly[2] - 2 * ctrl_poly[1] + ctrl_poly[0]) * (6 * (1.0 - t)) + \
           (ctrl_poly[3] - 2 * ctrl_poly[2] + ctrl_poly[1]) * (6 * t)

# Cr�er un vecteur tangent unitaire entre deux points
def create_tangent(pointA, pointB):
    return normalize(np.subtract(pointA, pointB))

# Normaliser un vecteur pour obtenir un vecteur unitaire
def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v

def get_control_points(beziers):
    control_points = []
    for bez in beziers:
        control_points.append(bez[0])
        control_points.append(bez[1])
        control_points.append(bez[2])
        control_points.append(bez[3])
    return control_points
        

