"""
Functions dealing with basic geometric problems.
"""
import numpy as np

def dot_product_from_points(ax1,ay1,ax2,ay2,bx1,by1,bx2,by2,):
    vector_a = np.array((ax2-ax1,ay2-ay1))
    vector_b = np.array((bx2-bx1,by2-by1))
    # print(vector_a,vector_b)

    dot_product = vector_a[0]*vector_b[0]+vector_a[1]*vector_b[1]

    return dot_product

def is_in_rect(a,b,d,m):
    """
    Provides a boolean value of whether point m is contained within
    rectangle abcd. All inputs should be arrays in order (x,y) defining
    the applicable vertices of abcd and the location of m.
    """
    m = np.array(m)
    mx = m[:,0]
    my = m[:,1]

    amab = dot_product_from_points(a[0],a[1],mx,my,
                                   a[0],a[1],b[0],b[1])
    abab = dot_product_from_points(a[0],a[1],b[0],b[1],
                                   a[0],a[1],b[0],b[1])
    amad = dot_product_from_points(a[0],a[1],mx,my,
                                   a[0],a[1],d[0],d[1])
    adad = dot_product_from_points(a[0],a[1],d[0],d[1],
                                   a[0],a[1],d[0],d[1])
    
    # print(amab)

    condition1 = 0<=amab
    condition2 = amab<=abab
    condition3 = 0<=amad
    condition4 = amad<=adad

    result = np.logical_and(np.logical_and(condition1,condition2), np.logical_and(condition3,condition4))
    
    # result = np.where(((0 <= amab) and (amab <= abab)) and ((0 <= amad) and (amad <= adad)), True, False)
    
    # result = ((0 <= amab) and (amab <= abab)) and ((0 <= amad) and (amad <= adad))

    return result

if __name__ == '__main__':
    a = (0,0)
    b = (0,5)
    d = (5,0)
    m = (0,4.9)

    result = is_in_rect(a,b,d,m)

    print(result)