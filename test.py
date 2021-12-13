import numpy as np

G = 1.

def acc(i, r, m):
    """ Return acceleration for index i"""
    sum = np.zeros_like(r[0])
    for j in range(len(m)):
        if j == i:
            continue
        else:
            sum += G * m[j] / np.linalg.norm(r[j] - r[i])**3 * (r[j] - r[i])
    
    return sum

def dT(i, r, v, m):
    """Return cinetic energy derived, but index i"""
    sum = 0
    for j in range(len(m)):
        if j == i:
            continue
        else:
            a = acc(j,r,m)
            sum += m[j] * np.dot(a,v[j])
    return sum

def dU(r,m):
    """Return all potential energy"""
    sum = 0
    n = len(m)
    for i in range(n):
        for j in range(i+1,n):
            sum += G * m[i] * m[j] / np.linalg.norm(r[i] - r[j])**2
    
    return sum

def get_v(i, r, v, m):
    dE = dT(i, r, v, m) + dU(r, m)
    accI = acc(i, r, m)